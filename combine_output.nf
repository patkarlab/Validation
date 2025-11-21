#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { trimming } from './npm1_mrd.nf'
include { trimming_aviti } from './npm1_mrd.nf'

process Combine_all {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy'
	
	input:
	val (Sample)

	output:
	tuple val(Sample), file("${Sample}_Final_merged.xlsx")

	script:
	"""
	cp $PWD/Final_Output/${Sample}/${Sample}.counts.bed ${Sample}.counts.bed
	cp $PWD/Final_Output/${Sample}/${Sample}.xlsx ${Sample}.xlsx

	if [ -f "$PWD/Final_Output/${Sample}/${Sample}_getitd/itds_collapsed-is-same_is-similar_is-close_is-same_trailing.tsv" ]; then
		cp "$PWD/Final_Output/${Sample}/${Sample}_getitd/itds_collapsed-is-same_is-similar_is-close_is-same_trailing.tsv" "${Sample}_getitd.tsv"
	else
		python ${params.getitd_gen} -o "${Sample}_getitd.tsv"
	fi

	python3 ${params.combine_all} --getitd_out ${Sample}_getitd.tsv --cov_out ${Sample}.counts.bed --merged_out ${Sample}.xlsx  --output ${Sample}_Final_merged.xlsx
	"""

}

process FLT3_ITD_EXT {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy'
	container "zhkddocker/flt3_itd_ext:v1.1"
	input:
		tuple val (Sample), file(paired_forward), file(paired_reverse)	
	output:
		tuple val (Sample), file("*FLT3_ITD.vcf"), file("*FLT3_ITD_summary.txt")
	script:
	"""
	perl /biosoft/FLT3_ITD_ext/FLT3_ITD_ext.pl -a 0 -f1 ${paired_forward} -f2 ${paired_reverse} -o ./ -g hg19 -n amplicon
	"""
}

process FLT3_ITD_EXT_FORMAT {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), file(ext_vcf), file(ext_summary)
	output:
		tuple val (Sample), file("${Sample}_ext.tsv")
	script:
	"""
	flt3_ext_format.py -s ${ext_summary} -v ${ext_vcf} -o ${Sample}_ext.tsv
	"""
}

workflow COMBINE {
	samples_ch = Channel.fromPath(params.input).splitCsv().flatten()
	Combine_all(samples_ch)
}

workflow FLT3_ITD {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }
	
	main:
	trimming(samples_ch)
	FLT3_ITD_EXT(trimming.out)
	FLT3_ITD_EXT_FORMAT(FLT3_ITD_EXT.out)
}

workflow FLT3_ITD_AVITI {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }
	
	main:
	trimming_aviti(samples_ch)
	FLT3_ITD_EXT(trimming_aviti.out)
	FLT3_ITD_EXT_FORMAT(FLT3_ITD_EXT.out)
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
