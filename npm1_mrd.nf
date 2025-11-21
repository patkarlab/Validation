#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""

bedfile = file("${params.bedfile}.bed", checkIfExists: true)
illumina_adapters = file("${params.adaptors}", checkIfExists: true )
nextera_adapters = file("${params.nextera_adapters}", checkIfExists: true )
genome_loc = file("${params.genome}", checkIfExists: true )
genome_dir = file("${genome_loc.parent}", checkIfExists: true)
genome_fasta = file("${genome_loc.name}")
filt3r_reference = file("${params.filt3r_ref}", checkIfExists: true)
getitd_path = file("${params.get_itd_path}")
minimap_getitd_reference = file("${params.genome_minimap_getitd}", checkIfExists: true)
annovar_path = file("${params.annovarLatest_path}")


process TRIM { 
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), file(read1), file(read2)
		file (illumina_adapters)
		file (nextera_adapters)
	output:
		tuple val (Sample), file("*1P.fq.gz"), file("*2P.fq.gz")
	script:
	"""	
	trimmomatic PE \
	${read1} ${read2} \
	-threads ${task.cpus} \
	-baseout ${Sample}.fq.gz \
	ILLUMINACLIP:${illumina_adapters}:2:30:10:2:keepBothReads \
	ILLUMINACLIP:${nextera_adapters}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40	
	sleep 5s
	"""
}

process MAPBAM {
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern : '*bam*'
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file (trim1), file (trim2)
		file (genome_dir)
		file (genome_fasta)
	output:
		tuple val (Sample), file ("${Sample}.bam"), file ("${Sample}.bam.bai")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:300" \
		-t ${task.cpus} ${genome_dir}/${genome_fasta} ${trim1} ${trim2} | samtools sort -@ ${task.cpus} -o ${Sample}.bam -

	samtools index ${Sample}.bam > ${Sample}.bam.bai
	"""
}

process FILT3R {
	tag "${Sample}"	
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*filt3r_json.csv'
	input:
		tuple val (Sample), file(trim1), file(trim2)
		file (filt3r_reference)
		file (annovar_path)
	output:
		tuple val (Sample), file("*_filt3r.vcf"), file("*_filt3r_json.csv"), file("*_filt3r_out.csv")

	script:
	"""
	filt3r -k 12 --ref ${filt3r_reference} --sequences ${trim1},${trim2} --nb-threads 64 --vcf --out ${Sample}_filt3r.json
	convert_json_to_csv.py ${Sample}_filt3r.json ${Sample}_filt3r_json.csv
	perl ${annovar_path}/convert2annovar.pl -format vcf4 ${Sample}_filt3r.vcf --outfile ${Sample}.filt3r.avinput --withzyg --includeinfo
	perl ${annovar_path}/table_annovar.pl ${Sample}.filt3r.avinput --out ${Sample}.filt3r --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${annovar_path}/humandb/ --xreffile ${annovar_path}/example/gene_fullxref.txt

	# Check if the multianno file is empty
	if [[ ! -s ${Sample}.filt3r.hg19_multianno.csv ]]; then
		touch ${Sample}.filt3r__final.csv
		touch ${Sample}_filt3r_json_filtered.csv
		touch ${Sample}_filt3r_out.csv
	else
		somaticseqoutput-format_filt3r.py ${Sample}.filt3r.hg19_multianno.csv ${Sample}.filt3r__final.csv
		filter_json.py ${Sample}_filt3r_json.csv ${Sample}_filt3r_json_filtered.csv
		merge_filt3r_csvs.py ${Sample}.filt3r__final.csv ${Sample}_filt3r_json_filtered.csv ${Sample}_filt3r_out.csv
	fi
	"""
}

process MINIMAP_GETITD {
	tag "${Sample}"
	label 'process_medium'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		tuple val(Sample), file(read1), file(read2)
		file (getitd_path)
		file (minimap_getitd_reference)
	output:
		tuple val (Sample), path ("${Sample}_getitd")
	script:
	"""
	minimap2 -ax sr -t ${task.cpus} ${minimap_getitd_reference} ${read1} ${read2} > ${Sample}.sam
	${params.samtools} view -b -h ${Sample}.sam -o ${Sample}.bam
	${params.samtools} sort ${Sample}.bam -o ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam
	${params.samtools} view ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${getitd_path}/getitd.py -reference ${getitd_path}/anno/amplicon.txt -anno ${getitd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}

process GETITD {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
			tuple val (Sample), file(bam), file(bamBai)
			file (getitd_path)
	output:
			path "*_getitd"
	script:
	"""
	${params.samtools} view ${bam} -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${getitd_path}/getitd.py -reference ${getitd_path}/anno/amplicon.txt -anno ${getitd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}

process FLT3_ITD_EXT {
	tag "${Sample}"
	input:
		tuple val (Sample), file(trim1), file(trim2)	
	output:
		tuple val (Sample), file("*vcf")
	script:
	"""
	perl /biosoft/FLT3_ITD_ext/FLT3_ITD_ext.pl -a 0 -f1 ${trim1} -f2 ${trim2} -o ./ -g hg19 -n amplicon

	"""
}


process COVERAGE {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), file(bam), file(bamBai)
		file (bedfile)
	output:
		tuple val (Sample), file ("${Sample}.counts.bed")
	script:
	"""
	bedtools bamtobed -i ${bam} > ${Sample}.bed
	bedtools coverage -counts -a ${bedfile} -b ${Sample}.bed > ${Sample}.counts.bed
	"""
}

process VARSCAN {
	tag "${Sample}"	
	input:
		tuple val (Sample), file(bam), file (bamBai)
		file (genome_dir)
		file (genome_fasta)
		file (bedfile)
	output:
		tuple val(Sample), file("${Sample}.varscan.vcf")
	script:
	"""
	${params.samtools} mpileup -d 1000000 -B -A -a -l ${bedfile} -f ${genome_dir}/${genome_fasta} ${bam} > ${Sample}.mpileup
	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-reads2 1 --min-avg-qual 5 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}

process ANNOVAR {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.varscan.csv'
	input:
		tuple val (Sample), file(varscanVcf)
		file (annovar_path)
	output:
		 tuple val (Sample), file ("*.hg19_multianno.csv"), file("*.varscan.csv")
	script:
	"""
	perl ${annovar_path}/convert2annovar.pl -format vcf4 ${varscanVcf} --outfile ${Sample}.varscan.avinput --withzyg --includeinfo
	perl ${annovar_path}/table_annovar.pl ${Sample}.varscan.avinput --out ${Sample}.varscan --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${annovar_path}/humandb/ --xreffile ${annovar_path}/example/gene_fullxref.txt

	if [ -s ${Sample}.varscan.hg19_multianno.csv ]; then
		somaticseqoutput-format_v2_varscan.py ${Sample}.varscan.hg19_multianno.csv ${Sample}.varscan.csv
	else
		touch ${Sample}.varscan.csv
	fi
	sleep 5s
	"""
}

process COMBINE_ALL {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.xlsx'
	input:
		tuple val (Sample),  file(filt3r_vcf), file(filt3r_json_csv), file(filt3r_out_csv), file (multianno), file(varscan_csv), path (getitd_dir), file(ext_vcf), file (coverage) 
	output:
		tuple val (Sample), file("${Sample}_ext.tsv"), file("${Sample}.xlsx")
	script:
	"""
	flt3_ext_format.py -v ${ext_vcf} -o ${Sample}_ext.tsv

	if [ -s ${getitd_dir}/itds_collapsed-is-same_is-similar_is-close_is-same_trailing.tsv ]; then
		GETITD_OUT="${getitd_dir}/itds_collapsed-is-same_is-similar_is-close_is-same_trailing.tsv"
	else
		touch itds_collapsed-is-same_is-similar_is-close_is-same_trailing.tsv
		GETITD_OUT="itds_collapsed-is-same_is-similar_is-close_is-same_trailing.tsv"
	fi
	combine_all.py -s ${Sample} -j ${filt3r_json_csv} -f ${filt3r_out_csv} -v ${varscan_csv} -g \$GETITD_OUT -x ${Sample}_ext.tsv -c ${coverage} -o ${Sample}.xlsx
	"""
}

process COMBINE_NPM1 {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.xlsx'
	input:
		tuple val (Sample), file (varscan_multianno), file(varscan_csv), file(coverage)
	output:
		tuple val (Sample), file("${Sample}.xlsx")	
	script:
	"""
	merge-tsv_npm1.py ${Sample}.xlsx ${varscan_csv} ${coverage}
	"""	
}


workflow NPM1_MRD {
	Channel.fromPath(params.input)
		.splitCsv(header:false)
		.map { row ->
			def sample = row[0].trim()
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz", checkIfExists: false)

			if (!r1 && !r2) {
				r1 = file("${params.sequences}/${sample}*_R1.fastq.gz", checkIfExists: false)
				r2 = file("${params.sequences}/${sample}*_R2.fastq.gz", checkIfExists: false)
			}
			tuple(sample, r1, r2)
		}
		.set { fastq_ch }
	main:
	TRIM(fastq_ch, illumina_adapters, nextera_adapters)
	MAPBAM(TRIM.out, genome_dir, genome_fasta)
	COVERAGE(MAPBAM.out, bedfile)
	VARSCAN(MAPBAM.out, genome_dir, genome_fasta, bedfile)
	ANNOVAR(VARSCAN.out, annovar_path)
	COMBINE_NPM1(ANNOVAR.out.join(COVERAGE.out))
}
workflow FLT3_MRD {
	Channel.fromPath(params.input)
		.splitCsv(header:false)
		.map { row ->
			def sample = row[0].trim()
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz", checkIfExists: false)

			if (!r1 && !r2) {
				r1 = file("${params.sequences}/${sample}*_R1.fastq.gz", checkIfExists: false)
				r2 = file("${params.sequences}/${sample}*_R2.fastq.gz", checkIfExists: false)
			}
			tuple(sample, r1, r2)
		}
		.set { fastq_ch }
	main:
	TRIM(fastq_ch, illumina_adapters, nextera_adapters)
	MAPBAM(TRIM.out, genome_dir, genome_fasta)
	FILT3R(TRIM.out, filt3r_reference, annovar_path)
	MINIMAP_GETITD(fastq_ch, getitd_path, minimap_getitd_reference)
	FLT3_ITD_EXT(TRIM.out)
	COVERAGE(MAPBAM.out, bedfile)
	VARSCAN(MAPBAM.out, genome_dir, genome_fasta, bedfile)
	ANNOVAR(VARSCAN.out, annovar_path)
	COMBINE_ALL(FILT3R.out.join(ANNOVAR.out.join(MINIMAP_GETITD.out.join(FLT3_ITD_EXT.out.join(COVERAGE.out)))))
}



workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
	def msg = """\
	Pipeline execution summary
	---------------------------
	Completed at : ${workflow.complete}
	Duration     : ${workflow.duration}
	""".stripIndent()

	println ""
	println msg
}
