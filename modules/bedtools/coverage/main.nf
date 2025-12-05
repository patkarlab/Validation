#!/usr/bin/nextflow

process COVERAGE {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), path (bam), path (bai)
		path (bedfile)
		val (bamtype)
	output:
		tuple val (Sample), path ("${Sample}_${bamtype}_coverage.bed")
	script:
	"""
	bedtools coverage -counts -a ${bedfile} -b ${bam} > "${Sample}_${bamtype}_coverage.bed"
	"""
}
