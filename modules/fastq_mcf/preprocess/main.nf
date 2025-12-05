#!/usr/bin/nextflow

process PREPROCESS {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), file(read1), file(read2)
		path(Adapt)
	output:
		tuple val(Sample), file("${Sample}_trim_R1.fastq.gz"), file("${Sample}_trim_R2.fastq.gz")
	script:
	"""	
	/opt/ea-utils/clipper/fastq-mcf -o ${Sample}_trim_R1.fastq.gz -o ${Sample}_trim_R2.fastq.gz -k 0 -q 20 ${Adapt} ${read1} ${read2}
	"""
}
