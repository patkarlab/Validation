#!/usr/bin/nextflow

process COMBINE_VCF {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), file(snp_vcf), file(indel_vcf)
	output:
		tuple val(Sample), file("${Sample}_varscan.vcf")
	script:
	""" 
	bzip2 -c ${snp_vcf} > ${snp_vcf}.gz
	bzip2 -c ${indel_vcf} > ${indel_vcf}.gz
	bcftools index -t ${snp_vcf}.gz
	bcftools index -t ${indel_vcf}.gz
	bcftools concat -a ${snp_vcf}.gz ${indel_vcf}.gz -o ${Sample}_varscan.vcf
	"""
}