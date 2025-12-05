#!/usr/bin/nextflow

process BAMTOFASTQ {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), path(bam)
	output:
		tuple val(Sample), path("${Sample}_R1.fastq.gz"), path("${Sample}_R2.fastq.gz")
	script:
	"""
	samtools fastq -@ ${task.cpus} -1 ${Sample}_R1.fastq.gz -2 ${Sample}_R2.fastq.gz -n ${bam}
	""" 
}