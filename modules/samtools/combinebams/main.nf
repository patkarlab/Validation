#!/usr/bin/nextflow

process COMBINEBAMS {
	tag "${Sample}"
	label 'process_low'
	input:
	        tuple val (Sample), path (cons_bam_list)
	output:
	        tuple val (Sample), path ("${Sample}.bam"), emit : combined_bam_ch
	script:
	"""
	samtools merge -@ ${task.cpus} ${Sample}.bam ${cons_bam_list}
	"""
}