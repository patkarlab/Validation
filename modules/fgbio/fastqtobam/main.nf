#!/usr/bin/nextflow

process ADD_UMI {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), path(read1), path (read2), path(umi)
	output:
		tuple val(Sample), path("${Sample}.unmapped.bam")
	script:
	"""
	fgbio -Xmx${task.memory.toGiga()}g --tmp-dir=. --async-io=true \
	FastqToBam --input ${read1} ${umi} ${read2} --read-structures +T +M +T --sample ${Sample} --library ${Sample} \
	--output ${Sample}.unmapped.bam
	"""
}
