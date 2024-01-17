#! /usr/bin/bash

samplesheet=$1
inputPath=$2
outputPath=$3

# ./creatingLink.sh samplesheet.csv /home/diagnostics/pipelines/Validation/sequences ./

for samples in `cat ${samplesheet}`
do
		ln -s "${inputPath}/${samples}_R1_001.fastq.gz" "${outputPath}"
		ln -s "${inputPath}/${samples}_R2_001.fastq.gz" "${outputPath}"
done
