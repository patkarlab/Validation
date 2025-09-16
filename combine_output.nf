#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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

workflow COMBINE {
	samples_ch = Channel.fromPath(params.input).splitCsv().flatten()
	Combine_all(samples_ch)
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
