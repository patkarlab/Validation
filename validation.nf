#!/usr/bin/env nextflow
nextflow.enable.dsl=2

"mkdir Coverview".execute()

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=

Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""
process trimming_trimmomatic {
	input:
		val Sample
	output:
		tuple val (Sample), file("*1P.fq.gz"), file("*2P.fq.gz")
	script:
	"""
	trimmomatic PE \
	${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz \
	-baseout ${Sample}.fq.gz \
	ILLUMINACLIP:${params.adaptors}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40
	sleep 5s
	"""
}

process minimap_getitd {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		val (Sample)
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

process mapping_reads{
	input:
		tuple val (Sample), file (pairAssembled)
	output:
		tuple val (Sample), file ("*.sam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" -M -t 20 ${params.genome} ${pairAssembled} > ${Sample}.sam
	"""
} 

process unpaird_mapping_reads{
	input:
		tuple val (Sample), file(paired_forward), file(paired_reverse)
	output:
		tuple val (Sample), file ("*.sam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" -M -t 20 ${params.genome} ${paired_forward} ${paired_reverse} > ${Sample}.sam
	"""
}

process sam_conversion{
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.sorted.bam*'
	input:
		tuple val (Sample), file (samfile)
	output:
		tuple val(Sample), file ("*.sorted.bam"), file ("*.sorted.bam.bai")
	script:
	"""
	${params.samtools} view -bT ${params.genome} ${samfile} > ${Sample}.bam
	${params.samtools} sort ${Sample}.bam > ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam > ${Sample}.sorted.bam.bai
	"""
}

process sam_conver_unpaired{
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.sorted.bam*'
	input:
		tuple val (Sample), file (samfile)
	output:
		tuple val(Sample), file ("*.sorted.bam"), file ("*.sorted.bam.bai")
	script:
	"""
	${params.samtools} view -bT ${params.genome} ${samfile} > ${Sample}.bam
	${params.samtools} sort ${Sample}.bam > ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam > ${Sample}.sorted.bam.bai
	"""
}		

process RealignerTargetCreator {
	input:
		tuple val (Sample), file (bamFile), file(bamBai)
	output:
		tuple val (Sample), file ("*.intervals")
	script:
	"""
	${params.java_path}/java -jar ${params.GATK38_path} -T RealignerTargetCreator -R ${params.genome} -nt 10 -I ${bamFile} --known ${params.site1} -o ${Sample}.intervals
	"""
}

process IndelRealigner{
	input:
		tuple val(Sample), file (targetIntervals), file(bamFile), file(bamBai)
	output:
		tuple val(Sample), file ("*.realigned.bam")
	script:
	"""
	echo ${Sample} ${targetIntervals} ${bamFile}
	${params.java_path}/java -jar ${params.GATK38_path} -T IndelRealigner -R ${params.genome} -I ${bamFile} -known ${params.site1} --targetIntervals ${targetIntervals} -o ${Sample}.realigned.bam
	"""
}

process BaseRecalibrator{
	input:
		tuple val (Sample), file (realignedBam)
	output:
		tuple val(Sample), file ("*.recal_data.table")
	script:
	"""
	${params.java_path}/java -jar ${params.GATK38_path} -T BaseRecalibrator -R ${params.genome} -I ${realignedBam} -knownSites ${params.site2} -knownSites ${params.site3} -maxCycle 600 -o ${Sample}.recal_data.table
	"""
}

process PrintReads{
	input:
		tuple val (Sample), file (realignedBam), file (recal_dataTable)
	output:
		tuple val (Sample), file ("*.aligned.recalibrated.bam")
	script:
	"""
	${params.java_path}/java -jar ${params.GATK38_path} -T PrintReads -R ${params.genome} -I ${realignedBam} --BQSR ${recal_dataTable} -o ${Sample}.aligned.recalibrated.bam
	"""
}

process generatefinalbam{
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam*'
	input:
		tuple val (Sample), file(alignedRecalibratedBam)
	output:
		tuple val(Sample), file ("*.final.bam"), file ("*.final.bam.bai"), file ("*.old_final.bam"), file ("*.old_final.bam.bai")
	script:
	"""
	${params.bedtools} sort -i ${params.bedfile}.bed > sorted.bed

	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${alignedRecalibratedBam} --out ${Sample}.abra.bam --ref ${params.genome} --threads 8 --targets sorted.bed --tmpdir ./ > abra.log

	${params.samtools} sort ${alignedRecalibratedBam} > ${Sample}.old_final.bam
	${params.samtools} index ${Sample}.old_final.bam > ${Sample}.old_final.bam.bai
	${params.samtools} sort ${Sample}.abra.bam > ${Sample}.final.bam
	${params.samtools} index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""
}

process hsmetrics_run{
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_hsmetrics.txt'
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*_hsmetrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${finalBam} O= ${Sample}_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	${params.hsmetrics_all} $PWD/Final_Output/hsmetrics.tsv ${Sample} ${Sample}_hsmetrics.txt
	"""
}

process InsertSizeMetrics {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*insert_size*' 
	input:
		tuple val(Sample), file (bamFile), file(bamBai)
	output:
		tuple val (Sample), file ("*_insert_size_metrics.txt"), file ("*_insert_size_histogram.pdf")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectInsertSizeMetrics I= ${bamFile} O= ${Sample}_insert_size_metrics.txt H= ${Sample}_insert_size_histogram.pdf M=0.5
	sleep 1s
	"""
}

process coverage {
	publishDir "$PWD/${Sample}/coverage/", mode: 'copy', pattern: '*.counts.bed'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*'
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("${Sample}.counts.bed"), file ("${Sample}_pindel.counts.bed")
	script:
	"""
	${params.bedtools} bamtobed -i ${finalBams[0]} > ${Sample}.bed
	${params.bedtools} coverage -counts -a ${params.bedfile}.bed -b ${Sample}.bed > ${Sample}.counts.bed
	${params.bedtools} coverage -counts -a ${params.flt3_bedfile}.bed -b ${Sample}.bed > ${Sample}_pindel.counts.bed
	mkdir -p $PWD/${Sample}/coverage/
	cp *.counts.bed $PWD/${Sample}/coverage/
	"""
}

process coverview_run {
	executor="local"
	publishDir "$PWD/${Sample}/Coverview/", mode: 'copy', pattern: '*.coverview_regions.csv'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.coverview_regions.csv'
	publishDir "$PWD/Coverview/", mode: 'copy', pattern: '*.coverview_regions.csv'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("${Sample}.coverview_regions.csv")
	script:
	"""
	${params.coverview_path}/coverview -i ${finalBam} -b ${params.bedfile}.bed -c ${params.coverview_path}/config/config.txt -o ${Sample}.coverview
	python3 ${params.coverview_script_path} ${Sample}.coverview_regions.txt ${Sample}.coverview_regions.csv
	"""
}

process coverview_report {
	executor="local"
	input:
		tuple val (Sample), file (coverview_csv)
	output:
		val Sample
	script:
	"""
	python3 ${params.coverview_report_path} ${PWD}/Coverview/ ${PWD}/Final_Output/
	"""
}

process mutect2_run{
	maxForks 10
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.mutect2.vcf'
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.mutect2.vcf")
	script:
	"""
	${params.samtools} view -bs 40.1 ${finalBam} > subsampled_01.bam
	${params.samtools} index subsampled_01.bam
	${params.mutect2} ${params.java_path} ${params.GATK38_path} ${params.genome} subsampled_01.bam ${Sample}.mutect2.vcf ${params.site2} ${params.bedfile}.bed	
	"""
}

process freebayes_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.freebayes.vcf'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.freebayes.vcf")
	script:
	"""
	${params.freebayes_path} -f ${params.genome} -b ${finalBam} -t ${params.bedfile}.bed > ${Sample}.freebayes.vcf
	"""
}

process vardict_run {
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.vardict.vcf'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	VarDict -G ${params.genome} -f 0.03 -N ${Sample} -b ${finalBam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.03 > ${Sample}.vardict.vcf	
	"""
}

process varscan_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan_snp.vcf'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan_indel.vcf'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan_snp.vcf.gz'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan_indel.vcf.gz'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan.vcf'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val(Sample), file("*.varscan.vcf")
	script:
	"""
	${params.samtools} mpileup -f ${params.genome} ${finalBam} > ${Sample}.mpileup
	${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}

process lofreq_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.lofreq.filtered.vcf'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val(Sample), file ("*.lofreq.filtered.vcf")
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample}.lofreq.pre.bam ${oldfinalBam}
	${params.samtools} sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	${params.lofreq_path} call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam
	${params.lofreq_path} filter -a 0.005 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	"""
}

process strelka_run{
	publishDir "$PWD/${Sample}/variants/strelka", mode: 'copy', pattern: '*.strelka.vcf'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple  val (Sample), file ("*.strelka.vcf")
	script:
	"""
	${params.strelka_path}/configureStrelkaGermlineWorkflow.py --bam ${finalBam} --referenceFasta ${params.genome} --callRegions  ${params.bedfile}.bed.gz --targeted --runDir ${PWD}/${Sample}/variants/strelka/
	${PWD}/${Sample}/variants/strelka/runWorkflow.py -m local -j 20
	gunzip -f ${PWD}/${Sample}/variants/strelka/results/variants/variants.vcf.gz
	mv ${PWD}/${Sample}/variants/strelka/results/variants/variants.vcf $PWD/${Sample}/variants/${Sample}.strelka.vcf
	cp $PWD/${Sample}/variants/${Sample}.strelka.vcf ./

	${params.strelka_path}/configureStrelkaSomaticWorkflow.py --normalBam ${params.NA12878_bam}  --tumorBam ${finalBam} --referenceFasta ${params.genome} --callRegions ${params.bedfile}.bed.gz --targeted --runDir ${PWD}/${Sample}/variants/strelka-somatic/
	${PWD}/${Sample}/variants/strelka-somatic/runWorkflow.py -m local -j 20
	
	${params.bcftools_path} concat -a ${PWD}/${Sample}/variants/strelka-somatic/results/variants/somatic.indels.vcf.gz ${PWD}/${Sample}/variants/strelka-somatic/results/variants/somatic.snvs.vcf.gz -o ${PWD}/${Sample}/variants/${Sample}.strelka-somatic.vcf
	"""
}

process platypus_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.platypus.vcf'
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val(Sample), file ("*.platypus.vcf")
	script:
	"""
	python2.7 ${params.platypus_path} callVariants --bamFiles=${finalBams[0]} --refFile=${params.genome} --output=${Sample}.platypus.vcf --nCPU=15 --minFlank=0 --filterDuplicates=0 --maxVariants=6 --minReads=6 --regions=${params.bedfile}_regions.txt
	"""
}

process somaticSeq_run {
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.somaticseq.vcf'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.avinput'
	publishDir "$PWD/${Sample}/ANNOVAR/", mode: 'copy', pattern: '*.hg19_multianno.csv'
	input:
			tuple val (Sample), file(freebayesVcf), file(platypusVcf), file(mutectVcf), file(vardictVcf), file(varscanVcf), file(lofreqVcf), file(strelkaVcf),file(finalBam), file(finalBamBai), file(oldfinalBam), file(oldfinalBamBai)
	output:
			tuple val (Sample), file ("*.somaticseq.vcf"), file("*.hg19_multianno.csv")
	script:
	"""
	${params.vcf_sorter_path} ${freebayesVcf} ${Sample}.freebayes.sorted.vcf
	${params.vcf_sorter_path} ${platypusVcf} ${Sample}.platypus.sorted.vcf

	python3 ${params.splitvcf_path} -infile ${Sample}.platypus.sorted.vcf -snv ${Sample}_platypus_cnvs.vcf -indel ${Sample}_platypus_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.freebayes.sorted.vcf -snv ${Sample}_freebayes_cnvs.vcf -indel ${Sample}_freebayes_indels.vcf

	${params.vcf_sorter_path} ${Sample}_platypus_cnvs.vcf ${Sample}_platypus_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_platypus_indels.vcf ${Sample}_platypus_indels_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_cnvs.vcf ${Sample}_freebayes_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_indels.vcf ${Sample}_freebayes_indels_sort.vcf

	somaticseq_parallel.py --output-directory ${PWD}/${Sample}/variants/${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --algorithm xgboost  --dbsnp-vcf  /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file $PWD/Final_Output/${Sample}/${Sample}.final.bam --mutect2-vcf ${PWD}/${Sample}/variants/${Sample}.mutect2.vcf --vardict-vcf ${PWD}/${Sample}/variants/${Sample}.vardict.vcf --varscan-vcf ${PWD}/${Sample}/variants/${Sample}.varscan.vcf --lofreq-vcf ${PWD}/${Sample}/variants/${Sample}.lofreq.filtered.vcf --strelka-vcf ${PWD}/${Sample}/variants/${Sample}.strelka.vcf --sample-name ${Sample} --arbitrary-snvs ${Sample}_freebayes_cnvs_sort.vcf ${Sample}_platypu_cnvs_sort.vcf --arbitrary-indels ${Sample}_freebayes_indels_sort.vcf ${Sample}_platypus_indels_sort.vcf

	${params.vcf_sorter_path} ${PWD}/${Sample}/variants/${Sample}.somaticseq/Consensus.sSNV.vcf ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf
	
	bgzip -c ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf > ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.vcf_sorter_path} ${PWD}/${Sample}/variants/${Sample}.somaticseq/Consensus.sINDEL.vcf ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf > ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf.gz

	${params.bcftools_path} concat -a ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf.gz ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.somaticseq.vcf

	sed -i 's/##INFO=<ID=MVDLK01,Number=7,Type=Integer,Description="Calling decision of the 7 algorithms: MuTect, VarScan2, VarDict, LoFreq, Strelka, SnvCaller_0, SnvCaller_1">/##INFO=<ID=MVDLKFP,Number=7,Type=String,Description="Calling decision of the 7 algorithms: MuTect, VarScan2, VarDict, LoFreq, Strelka, Freebayes, Platypus">/g' ${Sample}.somaticseq.vcf

	sed -i 's/MVDLK01/MVDLKFP/g' ${Sample}.somaticseq.vcf
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.somaticseq.vcf  --outfile ${Sample}.somaticseq.avinput --withzyg --includeinfo
	cp ${Sample}.somaticseq.vcf ${PWD}/Final_Output/${Sample}/
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.somaticseq.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	${params.cancervar} ${Sample}.somaticseq.hg19_multianno.csv ${Sample}
	cp ${Sample}myanno.hg19_multianno.txt.cancervar.ensemble.pred $PWD/${Sample}/
	"""
}

process mocha {
	publishDir "$PWD/Final_Output/${Sample}/MoChA/", mode: 'copy', pattern: '*.png'
	publishDir "$PWD/Final_Output/${Sample}/MoChA/", mode: 'copy', pattern: '*.pdf'
	publishDir "$PWD/Final_Output/${Sample}/MoChA/", mode: 'copy', pattern: '*.tsv'
	input:
		tuple val(Sample), file (somaticseqVcf), file (MultiannoCsv)
	output:
		tuple val(Sample), file ("*.png"), file ("*.pdf"), file ("*.tsv")
	script:
	"""
	mv ${somaticseqVcf} ${Sample}.somaticseq_old.vcf
	${params.bedtools} intersect -a ${Sample}.somaticseq_old.vcf -b ${params.mocha_bedfile} -header > ${Sample}.somaticseq.vcf
	${params.mocha} ${Sample} ./
	"""
}
process igv_reports{
	input:
		tuple val(Sample), file (somaticVcf), file (somaticseqMultianno)	
	output:
		val(Sample)
	script:
	"""
	perl ${params.annovarLatest_path}/table_annovar.pl ${somaticVcf} --out ${Sample}.annovar --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring . --otherinfo --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt -vcfinput

	${params.igv_script} ${params.genome} ${Sample}.annovar.hg19_multianno.vcf $PWD/Final_Output/${Sample}/${Sample}.final.bam $PWD/Final_Output/${Sample}/${Sample}_igv.html
	sleep 1s
	"""
}

process cnvkit_run {
	publishDir "$PWD/${Sample}/cnvkit/", mode: 'copy', pattern: '*.cn*'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		val Sample
	script:
	"""
	if [ ! -d ${PWD}/${Sample}/cnvkit ];then
		mkdir -p ${PWD}/${Sample}/cnvkit
	fi
	${params.cnvkit_path} ${finalBam} ${params.cnvkitRef} ${PWD}/${Sample}/cnvkit/
	#/${params.gene_scatter}/custom_scatter_v3.py ${params.gene_scatter}/chr_list_all.txt ${PWD}/${Sample}/cnvkit/${Sample}.final.cnr ${PWD}/${Sample}/cnvkit/${Sample}.final.cns ${Sample}
	/${params.gene_scatter}/custom_scatter_chrwise.py ${params.gene_scatter_list}/chrwise_list.txt ${PWD}/${Sample}/cnvkit/${Sample}.final.cnr ${PWD}/${Sample}/cnvkit/${Sample}.final.cns ${Sample}_chr_
	cp *gene_scatter.pdf $PWD/${Sample}/cnvkit/
	cp *gene_scatter.pdf $PWD/Final_Output/${Sample}/
	cp ${PWD}/${Sample}/cnvkit/${Sample}.final-scatter.png ${PWD}/${Sample}/cnvkit/${Sample}.final-diagram.pdf ${PWD}/Final_Output/${Sample}/	
	"""
}

process ifcnv_run {
	//publishDir "$PWD/Final_Output/ifCNV/", mode: 'copy', pattern: '*.html'
	//publishDir "$PWD/Final_Output/ifCNV/", mode: 'copy', pattern: '*.tsv'
	input:
		val Sample
	output:
		val (Sample)
	script:
	"""
	${params.links} $PWD/Final_Output/ ${params.input}
	#${params.ifcnv} ./ ${params.bedfile}.bed ./
	#touch -c *.bam.bai
	mkdir ifCNV_chrwise
	mkdir ifCNV_genewise
	${params.ifcnv} ./ ${PWD}/bedfiles/CNV_small_hg19renamed_sortd.bed ifCNV_chrwise
	${params.ifcnv} ./ ${PWD}/bedfiles/newmyeloid_sortd.bed ifCNV_genewise

	for i in `cat ${params.input}`
	do
		if [ ! -d $PWD/Final_Output/\${i}/ifCNV_chrwise ]; then
			mkdir $PWD/Final_Output/\${i}/ifCNV_chrwise
		fi

		if [ ! -d $PWD/Final_Output/\${i}/ifCNV_genewise ]; then
			mkdir $PWD/Final_Output/\${i}/ifCNV_genewise
		fi
	done

	if [ -f ifCNV_chrwise/ifCNV.tsv ]; then
		for i in `awk 'NR>1{print \$3}' ifCNV_chrwise/ifCNV.tsv | awk 'BEGIN{FS="."}{print \$1}' | sort | uniq`
		do
			cp ifCNV_chrwise/\${i}.*.html $PWD/Final_Output/\${i}/ifCNV_chrwise/
		done
	fi

	if [ -f ifCNV_genewise/ifCNV.tsv ]; then
		for i in `awk 'NR>1{print \$3}' ifCNV_genewise/ifCNV.tsv | awk 'BEGIN{FS="."}{print \$1}' | sort | uniq` 
		do		
			cp ifCNV_genewise/\${i}.*.html $PWD/Final_Output/\${i}/ifCNV_genewise/
		done
	fi
	"""
}

process ifcnv {
	//publishDir "$PWD/Final_Output/ifCNV/", mode: 'copy', pattern: '*.html'
	//publishDir "$PWD/Final_Output/ifCNV/", mode: 'copy', pattern: '*.tsv'
	//input:
		//val Sample
	//output:
		//tuple val (Sample)
	script:
	"""
	${params.links} $PWD/sequences/ ${params.input}
	#${params.ifcnv} ./ ${params.bedfile}.bed ./
	#touch -c *.bam.bai
	mkdir ifCNV_chrwise
	mkdir ifCNV_genewise
	${params.ifcnv} ./ ${PWD}/bedfiles/CNV_small_hg19renamed_sortd.bed ifCNV_chrwise
	${params.ifcnv} ./ ${PWD}/bedfiles/newmyeloid_sortd.bed ifCNV_genewise

	for i in `cat ${params.input}`
	do
		if [ ! -d $PWD/Final_Output/\${i}/ifCNV_chrwise ]; then
			mkdir $PWD/Final_Output/\${i}/ifCNV_chrwise
		fi

		if [ ! -d $PWD/Final_Output/\${i}/ifCNV_genewise ]; then
			mkdir $PWD/Final_Output/\${i}/ifCNV_genewise
		fi
	done

	if [ -f ifCNV_chrwise/ifCNV.tsv ]; then
		for i in `awk 'NR>1{print \$3}' ifCNV_chrwise/ifCNV.tsv | awk 'BEGIN{FS="."}{print \$1}' | sort | uniq`
		do
			cp ifCNV_chrwise/\${i}.*.html $PWD/Final_Output/\${i}/ifCNV_chrwise/
		done
	fi

	if [ -f ifCNV_genewise/ifCNV.tsv ]; then
		for i in `awk 'NR>1{print \$3}' ifCNV_genewise/ifCNV.tsv | awk 'BEGIN{FS="."}{print \$1}' | sort | uniq`
		do
			cp ifCNV_genewise/\${i}.*.html $PWD/Final_Output/\${i}/ifCNV_genewise/
		done
	fi
	"""
}

process update_db {
	publishDir "$PWD/work/", mode: 'copy', pattern: 'freq.txt'
	input:
		val Sample
	output:
		file ("freq.txt")
	script:
	"""
	for i in `cat ${params.input}`
	do 
		if [ -f ${PWD}/Final_Output/\${i}/\${i}.somaticseq.vcf ]; then
			ln -s ${PWD}/Final_Output/\${i}/\${i}.somaticseq.vcf ./
		fi
	done
	files=\$(ls *.somaticseq.vcf)
	${params.updatedb} ${params.icmrdb} \${files}
	${params.freq_py} ${params.icmrdb} freq.txt
	"""
}

process pindel {
	publishDir "$PWD/${Sample}/pindel/", mode: 'copy', pattern: '*pindel_SI.vcf'
	publishDir "$PWD/${Sample}/pindel/", mode: 'copy', pattern: '*.avinput'
	publishDir "$PWD/${Sample}/pindel/", mode: 'copy', pattern: '*_pindel.hg19_multianno.csv'
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file("*pindel.hg19_multianno.csv"), file("*.avinput"), file("*pindel_SI.vcf")
	script:
	"""
	export BAM_2_PINDEL_ADAPT=${params.pindel}/Adaptor.pm
	sh ${params.pindel_config_script} -s ${Sample} -b ${finalBam} -c config.txt
	${params.pindel}/pindel -f ${params.genome} -i config.txt -c chr13 -o ${Sample}_pindel
	${params.pindel}/pindel2vcf -r ${params.genome} -P ${Sample}_pindel -R hg19 -d 07102019 -v ${Sample}_pindel_SI.vcf

	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}_pindel_SI.vcf --outfile ${Sample}_pindel.avinput --withzyg --includeinfo

	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}_pindel.avinput ${params.annovarLatest_path}/humandb/ -buildver hg19 -out ${Sample}_pindel --remove -protocol refGene,cytoBand,cosmic84 --operation g,r,f -nastring '.' --otherinfo --csvout --thread 10 --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	"""
}

process combine_variants{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.avinput'
	publishDir "$PWD/${Sample}/ANNOVAR/", mode: 'copy', pattern: '*.hg19_multianno.csv'
	input:
		tuple val (Sample), file(freebayesVcf), file(platypusVcf)
	output:
		tuple val(Sample), file ("*.combined.vcf"),  file ("*.hg19_multianno.csv")
	script:
	"""
	grep "^#" ${PWD}/${Sample}/variants/${Sample}.freebayes.vcf > ${Sample}.freebayes.sorted.vcf
	grep -v "^#" ${PWD}/${Sample}/variants/${Sample}.freebayes.vcf | sort -k1,1V -k2,2g >> ${Sample}.freebayes.sorted.vcf
	grep "^#" ${PWD}/${Sample}/variants/${Sample}.platypus.vcf > ${Sample}.platypus.sorted.vcf
	grep -v "^#" ${PWD}/${Sample}/variants/${Sample}.platypus.vcf | sort -k1,1V -k2,2g >> ${Sample}.platypus.sorted.vcf
	${params.java_path}/java -jar ${params.GATK38_path} -T CombineVariants -R ${params.genome} --variant ${Sample}.freebayes.sorted.vcf --variant ${Sample}.platypus.sorted.vcf -o ${Sample}.combined.vcf -genotypeMergeOptions UNIQUIFY
	
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.combined.vcf  --outfile ${Sample}.combined.avinput --withzyg --includeinfo

	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.combined.avinput --out ${Sample}.combined --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	"""
}

process cava {
	publishDir "$PWD/${Sample}/CAVA/", mode: 'copy'
	input:
		tuple val(Sample), file (somaticVcf), file (somaticseqMultianno)
	output:
		tuple val(Sample), file ("*.cava.csv")
	script:
	"""
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i ${somaticVcf} -o ${Sample}.somaticseq
	python3 ${params.cava_script_path} ${Sample}.somaticseq.txt ${Sample}.cava.csv
	"""
}

process format_somaticseq_combined {
	input:
		tuple val (Sample), file(somaticseqVcf), file (multianno), file (combinedVcf),  file (hg19_multianno)
	output:
		val Sample
	script:
	"""
	mkdir "$PWD/${Sample}/Annovar_Modified/"
	python3 ${params.format_somaticseq_script} ${PWD}/${Sample}/ANNOVAR/${Sample}.somaticseq.hg19_multianno.csv ${PWD}/${Sample}/Annovar_Modified/${Sample}.somaticseq.csv
	python3 ${params.format_combined_script} ${PWD}/${Sample}/ANNOVAR/${Sample}.combined.hg19_multianno.csv ${PWD}/${Sample}/Annovar_Modified/${Sample}.combined.csv
	"""
}

process format_concat_combine_somaticseq {
	input:
		val (Sample)
	output:
		val Sample
	script:
	"""
	sed -i '1d' ${PWD}/${Sample}/Annovar_Modified/${Sample}.combined.csv
	sed -i '1d' ${PWD}/${Sample}/Annovar_Modified/${Sample}.somaticseq.csv

	cp ${PWD}/${Sample}/Annovar_Modified/${Sample}.somaticseq.csv ${PWD}/${Sample}/Annovar_Modified/${Sample}.concat.csv
	python3 ${params.format_remove_artefact_script} ${PWD}/${Sample}/Annovar_Modified/${Sample}.concat.csv ${params.artefactFile} ${PWD}/${Sample}/Annovar_Modified/${Sample}.final.concat.csv ${PWD}/${Sample}/Annovar_Modified/${Sample}.artefacts.csv

	sed -i '1iChr,Start,End,Ref,Alt,Variant_Callers,FILTER,SOMATIC_FLAG,VariantCaller_Count,REF_COUNT,ALT_COUNT,VAF,Func.refGene,Gene.refGene,ExonicFunc.refGene,AAChange.refGene,Gene_full_name.refGene,Function_description.refGene,Disease_description.refGene,cosmic84,PopFreqMax,1000G_ALL,ExAC_ALL,CG46,ESP6500siv2_ALL,InterVar_automated' ${PWD}/${Sample}/Annovar_Modified/${Sample}.final.concat.csv
	
	sed -i '1iChr,Start,End,Ref,Alt,Variant_Callers,FILTER,SOMATIC_FLAG,VariantCaller_Count,REF_COUNT,ALT_COUNT,VAF,Func.refGene,Gene.refGene,ExonicFunc.refGene,AAChange.refGene,Gene_full_name.refGene,Function_description.refGene,Disease_description.refGene,cosmic84,PopFreqMax,1000G_ALL,ExAC_ALL,CG46,ESP6500siv2_ALL,InterVar_automated' ${PWD}/${Sample}/Annovar_Modified/${Sample}.artefacts.csv
	"""
}

process format_pindel {
	//errorStrategy 'ignore'
	publishDir "$PWD/${Sample}/pindel/", mode: 'copy', pattern: '*_final.pindel.csv'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_final.pindel.csv'
	input:
		tuple val (Sample), file(pindel_hg19_multianno), file(avinput), file(pindel_SI_vcf), file (sample_counts_bed), file (sample_pindel_counts_bed)	
	output:
		val Sample
	script:
	"""
	python3 ${params.format_pindel_script} ${PWD}/${Sample}/coverage/${Sample}_pindel.counts.bed ${PWD}/${Sample}/pindel/${Sample}_pindel.hg19_multianno.csv ${PWD}/${Sample}/pindel/${Sample}_final.pindel.csv

	"""
}

process merge_csv {
	input:
		tuple val (Sample), file (cava_csv)
	output:
		tuple val (Sample), file ("output_temp.xlsx")
	script:
	"""
	sed -i 's/\t/,/g' ${PWD}/${Sample}/cnvkit/${Sample}.final.cnr
	if [ -f ${PWD}/${Sample}/Annovar_Modified/${Sample}.final.concat.csv ]; then
		python3 ${params.pharma_marker_script} ${Sample} ${PWD}/${Sample}/Annovar_Modified/ ${params.pharma_input_xlxs} ${PWD}/${Sample}/${Sample}_pharma.csv
	else
		touch ${PWD}/${Sample}/${Sample}_pharma.csv
	fi
	python3 ${params.merge_csvs_script} ${Sample} ${PWD}/${Sample}/Annovar_Modified/ ${PWD}/Final_Output/${Sample}/${Sample}.xlsx ${PWD}/${Sample}/CAVA/ ${PWD}/${Sample}/Coverview/${Sample}.coverview_regions.csv ${PWD}/${Sample}/pindel/${Sample}_final.pindel.csv ${PWD}/${Sample}/cnvkit/${Sample}.final.cnr ${PWD}/${Sample}/${Sample}_pharma.csv

	cp ${PWD}/${Sample}/Annovar_Modified/${Sample}.final.concat.csv ${Sample}.final.concat_append.csv
	${params.vep_script_path} ${PWD}/Final_Output/${Sample}/${Sample}.somaticseq.vcf ${PWD}/Final_Output/${Sample}/${Sample}
	${params.vep_extract_path} ${Sample}.final.concat_append.csv ${PWD}/Final_Output/${Sample}/${Sample}_vep_delheaders.txt > ${Sample}.vep
	${params.cancervar_extract} $PWD/${Sample}/${Sample}myanno.hg19_multianno.txt.cancervar.ensemble.pred ${Sample}.vep ${Sample}_cancervar.csv

	${params.pcgr_cpsr_script_path} ${PWD}/Final_Output/${Sample}/${Sample}.xlsx ${Sample}_cancervar.csv
	cp output_temp.xlsx ${PWD}/Final_Output/${Sample}/${Sample}.xlsx
	
	## for adding CNV_coverview_region and myeloid_coverview_region sheet in final xlsx  
	python3 ${params.myeloid_cnv} ${PWD}/Final_Output/${Sample}/${Sample}.coverview_regions.csv /home/pipelines/MMpanel/scripts/cnvkit_cnvmyeloid/CNV_3072_hg19_RegionNames.txt ${PWD}/Final_Output/${Sample}/${Sample}.xlsx
	"""
}

process update_freq {
	input:
		val (Sample)
	output:
		val (Sample)
	script:
	"""
	ln -s ${PWD}/work/freq.txt ./
	for i in `cat ${params.input}`
	do
		if [ -f ${PWD}/Final_Output/\${i}/\${i}.xlsx ]; then
			${params.update_freq_excel} ${PWD}/Final_Output/\${i}/\${i}.xlsx freq.txt
		fi
	done
	"""
}

process Final_Output {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.png'
	input:
		tuple val (Sample), file (sample_counts_bed), file (sample_pindel_counts_bed)
	output:
		tuple val (Sample), file ("*.png")
	script:
	"""
	python3 ${params.coveragePlot_script} ${Sample} $PWD/${Sample}/coverage/${Sample}.counts.bed ./
	"""
}
process remove_files{
	input:
		tuple val (Sample),file (output_excel), file (coverage_png)
	script:
	"""
	rm -rf ${PWD}/${Sample}/
	"""
}

workflow VALIDATION {
	
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }
	
	main:
	trimming_trimmomatic(samples_ch) | pair_assembly_pear | mapping_reads | sam_conversion
	//unpaird_mapping_reads(trimming_trimmomatic.out) | sam_conver_unpaired
	//minimap_getitd(samples_ch)
	RealignerTargetCreator(sam_conversion.out)
	IndelRealigner(RealignerTargetCreator.out.join(sam_conversion.out)) | BaseRecalibrator
	PrintReads(IndelRealigner.out.join(BaseRecalibrator.out)) | generatefinalbam
	hsmetrics_run(generatefinalbam.out)
	//InsertSizeMetrics(sam_conver_unpaired.out)
	//coverage(generatefinalbam.out)
	coverview_run(generatefinalbam.out)
	//coverview_report(coverview_run.out)
	//platypus_run(generatefinalbam.out)
	//freebayes_run(generatefinalbam.out)
	//mutect2_run(generatefinalbam.out)
	//vardict_run(generatefinalbam.out)
	//varscan_run(generatefinalbam.out)
	//lofreq_run(generatefinalbam.out)
	//strelka_run(generatefinalbam.out)
	//cnvkit_run(generatefinalbam.out)
	//pindel(generatefinalbam.out)
	//format_pindel(pindel.out.join(coverage.out))
	//somaticSeq_run(freebayes_run.out.join(platypus_run.out.join(mutect2_run.out.join(vardict_run.out.join(varscan_run.out.join(lofreq_run.out.join(strelka_run.out.join(generatefinalbam.out))))))))
	//igv_reports(somaticSeq_run.out)
	//ifcnv_run(generatefinalbam.out.collect())
	//update_db(somaticSeq_run.out.collect())
	//mocha(somaticSeq_run.out)
	//combine_variants(freebayes_run.out.join(platypus_run.out))
	//cava(somaticSeq_run.out)
	//format_somaticseq_combined(somaticSeq_run.out.join(combine_variants.out))
	//format_concat_combine_somaticseq(format_somaticseq_combined.out)	
	//merge_csv(format_concat_combine_somaticseq.out.join(cava.out.join(format_pindel.out.join(cnvkit_run.out))))
	//update_freq(merge_csv.out.collect())
	//Final_Output(coverage.out.join(cnvkit_run.out))
	//remove_files(merge_csv.out.join(Final_Output.out))
}

workflow IfCnv {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }
	
	main:
	ifcnv()
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
