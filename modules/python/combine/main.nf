#!/usr/bin/env nextflow

process COMBINE_CALLERS {
	tag "${Sample}"
	label 'process_single'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), path(mutect2), path(vardict), path(varscan), path(coverage)
		val (bamtype)
	output:
		tuple val (Sample), path ("${Sample}_${bamtype}.xlsx")
	script:
	"""
	somaticseqoutput-format_v2_mutect2.py ${mutect2} mutect2_.csv
	somaticseqoutput-format_v2_varscan.py ${varscan} varscan_.csv
	somaticseqoutput-format_v2_vardict.py ${vardict} vardict_.csv
	combine_callers.py ${Sample}_${bamtype}.xlsx mutect2_.csv varscan_.csv vardict_.csv ${coverage}	
	"""
}