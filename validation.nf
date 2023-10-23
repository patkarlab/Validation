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
	//publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.sorted.bam*'
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
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T RealignerTargetCreator -R ${params.genome} -nt 10 -I ${bamFile} --known ${params.site1} -o ${Sample}.intervals
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
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T IndelRealigner -R ${params.genome} -I ${bamFile} -known ${params.site1} --targetIntervals ${targetIntervals} -o ${Sample}.realigned.bam
	"""
}

process BaseRecalibrator{
	input:
		tuple val (Sample), file (realignedBam)
	output:
		tuple val(Sample), file ("*.recal_data.table")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T BaseRecalibrator -R ${params.genome} -I ${realignedBam} -knownSites ${params.site2} -knownSites ${params.site3} -maxCycle 600 -o ${Sample}.recal_data.table
	"""
}

process PrintReads{
	input:
		tuple val (Sample), file (realignedBam), file (recal_dataTable)
	output:
		tuple val (Sample), file ("*.aligned.recalibrated.bam")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T PrintReads -R ${params.genome} -I ${realignedBam} --BQSR ${recal_dataTable} -o ${Sample}.aligned.recalibrated.bam
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
	errorStrategy 'ignore'
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

process remove_files{
	errorStrategy 'ignore'
	input:
		tuple val (Sample),file (InsertSizeMetrics),file (InsertSizeMetricsPdf), file (hsmetrics), file (coverageBed), file (pindelcoveragebed), file (coverview)
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
	unpaird_mapping_reads(trimming_trimmomatic.out) | sam_conver_unpaired
	RealignerTargetCreator(sam_conversion.out)
	IndelRealigner(RealignerTargetCreator.out.join(sam_conversion.out)) | BaseRecalibrator
	PrintReads(IndelRealigner.out.join(BaseRecalibrator.out)) | generatefinalbam
	hsmetrics_run(generatefinalbam.out)
	InsertSizeMetrics(sam_conver_unpaired.out)
	coverage(generatefinalbam.out)
	coverview_run(generatefinalbam.out)
	//coverview_report(coverview_run.out)
	remove_files(InsertSizeMetrics.out.join(hsmetrics_run.out.join(coverage.out.join(coverview_run.out))))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
