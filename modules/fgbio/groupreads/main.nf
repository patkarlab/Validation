#!/usr/bin/nextflow

process GROUPREADSBYUMI {
	tag "${Sample}"
	maxForks 5
	label 'process_medium'
	// publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_grouped.bam'
	input:
		tuple val (Sample), path(mapped_bam)
	output:
		tuple val (Sample), path ("*_grouped.bam"), emit : grouped_bam_ch
	script:
	"""
	for bam_file in ${mapped_bam}
	do
		file_name=\$( basename \${bam_file} .bam)       # Removing the .bam extension
		fgbio -Xmx${task.memory.toGiga()}g --compression 1 --async-io GroupReadsByUmi --input \${bam_file} --strategy Adjacency --edits 1 \
		--output \${file_name}_grouped.bam --family-size-histogram \${file_name}.tag-family-sizes.txt
	done
	"""
}