#!/usr/bin/nextflow

process SPLITBAM {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), path (mapbam), path(bam_bai)
		path (bedfile)
	output:
		tuple val (Sample), path ("${Sample}_*_primary.bam"), emit: chrwise_bamlist
	script:
	"""
	samtools sort -@ ${task.cpus} ${mapbam} -o ${Sample}_sortd.bam
	samtools index ${Sample}_sortd.bam > ${Sample}_sortd.bam.bai
	samtools idxstats ${Sample}_sortd.bam | cut -f 1 | grep -v '*' > chromosomes.txt
	for i in `cat chromosomes.txt`
	do
		# Obtain the read names (headers) aligned to each chromosome
		samtools view -@ ${task.cpus} ${Sample}_sortd.bam "\${i}" | cut -f1 | sort -u > \${i}_names.txt
		# Mate rescue: extract *all* alignments for those names
		samtools view -@ ${task.cpus} -N \${i}_names.txt -b ${Sample}_sortd.bam > ${Sample}_\${i}_rescued.bam
		# 3) Keep only primaries (drop secondary=256 and supplementary=2048)
		samtools view -@ ${task.cpus} -F 2304 -b ${Sample}_\${i}_rescued.bam > ${Sample}_\${i}_primary.bam
	done
	"""	
}