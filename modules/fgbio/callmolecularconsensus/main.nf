#!/usr/bin/nextflow

process CALLMOLCONSREADS {
	tag "${Sample}"
	maxForks 5
	label 'process_high'
	input:
		tuple val (Sample), path (grouped_bam)
		path (GenFile)
		path (GenDir)
	output:
		tuple val (Sample), path ("*_cons_umap.bam"), emit : consensus_bam_ch
	script:
	"""
	# This step generates unmapped consensus reads from the grouped reads and immediately filters them
	for bams in ${grouped_bam}
	do
		outfile_name=\$( basename \${bams} .bam)        # Removing the .bam extension
		fgbio -Xmx${task.memory.toGiga()}g --compression 0 CallMolecularConsensusReads --input \${bams} \
		--output /dev/stdout --min-reads 2 --min-input-base-quality 20 --threads ${task.cpus} | fgbio -Xmx${task.memory.toGiga()}g --compression 1 \
		FilterConsensusReads --input /dev/stdin --output \${outfile_name}_cons_umap.bam --ref ${GenFile} \
		--min-reads 2 --min-base-quality 20
	done
	"""
}