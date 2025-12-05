#!/usr/bin/nextflow

process HSMETRICS {
	tag "${Sample}"
	label 'process_low'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val(Sample), path(bam), path(bai)
		path (bedfile)
		path (GenFile)
		path (GenDir)
		val (bamtype)
	output:
		tuple val(Sample), path("${Sample}_${bamtype}_hsmetrics.txt")
	script:
	"""
	gatk BedToIntervalList I=${bedfile} O=${bedfile}_sortd.interval_list SD=${GenFile}.dict
	gatk CollectHsMetrics I=${bam} O=${Sample}_${bamtype}_hsmetrics.txt BAIT_INTERVALS=${bedfile}_sortd.interval_list TARGET_INTERVALS=${bedfile}_sortd.interval_list R= ${GenFile} VALIDATION_STRINGENCY=LENIENT
	"""
}