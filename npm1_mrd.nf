#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""

process minimap_getitd {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
			val Sample
	output:
			path "*_getitd"
	script:
	"""
	minimap2 -ax sr ${params.genome_minimap_getitd} ${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz > ${Sample}.sam
	${params.samtools} view -b -h ${Sample}.sam -o ${Sample}.bam
	${params.samtools} sort ${Sample}.bam -o ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam
	${params.samtools} view ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${params.get_itd_path}/getitd.py -reference ${params.get_itd_path}/anno/amplicon.txt -anno ${params.get_itd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}

process trimming { 
	input:
		val (Sample)
	output:
		//tuple val (Sample), file("${Sample}.R1.trimmed.fastq"), file("${Sample}.R2.trimmed.fastq")
		tuple val (Sample), file("*1P.fq.gz"), file("*2P.fq.gz")
	script:
	"""	
	#/home/programs/ea-utils/clipper/fastq-mcf -o ${Sample}.R1.trimmed.fastq -o ${Sample}.R2.trimmed.fastq -l 53 -k 0 -q 0 /home/diagnostics/pipelines/smMIPS_pipeline/code/functions/preprocess_reads_miseq/smmip_adaptors.fa ${params.sequences}/${Sample}_S*_R1_*.fastq.gz  ${params.sequences}/${Sample}_S*_R2_*.fastq.gz
	trimmomatic PE \
	${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz \
	-baseout ${Sample}.fq.gz \
	ILLUMINACLIP:${params.adaptors}:2:30:10:2:keepBothReads \
	ILLUMINACLIP:${params.nextera_adapters}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40	
	sleep 5s
	"""
}

process pair_assembly_pear {
	input:
		tuple val (Sample), file(paired_forward), file(paired_reverse)
	output:
		tuple val (Sample), file("*.assembled.fastq")
	script:
	"""
	${params.pear_path} -f ${paired_forward} -r ${paired_reverse} -o ${Sample} -n 53 -j 25
	"""
}

process mapping_reads {
	input:
		tuple val (Sample), file (pairAssembled)
	output:
		tuple val (Sample), file ("*.sam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" -M -t 20 ${params.genome} ${pairAssembled} > ${Sample}.sam
	"""
}

process map_paired_reads {
	input:
		tuple val (Sample), file (paired_read1), file (paired_read2)
	output:
		tuple val (Sample), file ("*.sam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:300" -t 20 ${params.genome} ${paired_read1} ${paired_read2} > ${Sample}.sam
	"""
}

process sam_conversion {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.uncollaps.bam*'
	input:
		tuple val (Sample), file (samfile)	
	output:
		tuple val(Sample), file ("*.uncollaps.bam"), file ("*.uncollaps.bam.bai")
	script:
	"""
	${params.samtools} view -bT ${params.genome} ${samfile} > ${Sample}.bam
	${params.samtools} sort ${Sample}.bam > ${Sample}.uncollaps.bam
	${params.samtools} index ${Sample}.uncollaps.bam > ${Sample}.uncollaps.bam.bai
	"""
}

process coverage {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai)
	output:
		tuple val (Sample), file ("${Sample}.counts.bed")
	script:
	"""
	${params.bedtools} bamtobed -i ${finalBams[0]} > ${Sample}.bed
	${params.bedtools} coverage -counts -a ${params.bedfile}.bed -b ${Sample}.bed > ${Sample}.counts.bed
	"""
}

process varscan {
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val(Sample), file("${Sample}.varscan.vcf")
	script:
	"""
	${params.samtools} mpileup -d 1000000 -B -A -a -l ${params.bedfile}.bed -f ${params.genome} ${finalBam} > ${Sample}.mpileup
	${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-reads2 1 --min-avg-qual 5 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup  --min-reads2 1 --min-avg-qual 5 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}

process Annovar {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.varscan_.csv'
	input:
		tuple val (Sample), file(varscanVcf)
	output:
		 tuple val (Sample), file ("*.hg19_multianno.csv"), file("*.varscan_.csv")
	script:
	"""
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${varscanVcf}  --outfile ${Sample}.varscan.avinput --withzyg --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.varscan.avinput --out ${Sample}.varscan --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	python3 ${PWD}/scripts/somaticseqoutput-format_v2_varscan.py ${Sample}.varscan.hg19_multianno.csv ${Sample}.varscan_.csv
	#cp ${Sample}.varscan_.csv $PWD/Final_Output/${Sample}/
	sleep 5s
	"""
}

workflow NPM1 {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }
	
	main:
	//trimming(samples_ch) | pair_assembly_pear | mapping_reads | sam_conversion
	trimming(samples_ch) | map_paired_reads | sam_conversion
	minimap_getitd(samples_ch)
	coverage(sam_conversion.out)
	varscan(sam_conversion.out)
	Annovar(varscan.out)
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}