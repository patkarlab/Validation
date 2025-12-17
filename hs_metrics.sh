#! /usr/bin/bash

for i in `head -n20 /home/diagnostics/pipelines/Validation/samplesheet_superset.csv` 
do 
	#java -jar /home/programs/picard/build/libs/picard.jar CollectHsMetrics I=/home/diagnostics/pipelines/Validation/Final_Output/$i/$i".sorted.bam" O=/home/diagnostics/pipelines/Validation/Final_Output/$i/$i".hsmetrics.txt" R=/home/reference_genomes/hg38_broad/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta BAIT_INTERVALS=/home/diagnostics/pipelines/Validation/bedfiles/ALLPanel_Bed_Coordinates_hg38_140524_format_sortd.interval_list TARGET_INTERVALS=/home/diagnostics/pipelines/Validation/bedfiles/ALLPanel_Bed_Coordinates_hg38_140524_format_sortd.interval_list VALIDATION_STRINGENCY=LENIENT

	echo -ne $i'\t'; grep -v '#' /home/diagnostics/pipelines/Validation/MRD_IDT/$i/$i"_uncollaps_hsmetrics.txt"  | awk 'BEGIN{FS="\t"; OFS="\t"}NR==3{ print $7,$8}'

	#${params.coverview_path}/coverview -i ${finalBam} -b ${params.bedfile}.bed -c ${params.coverview_path}/config/config.txt -o ${Sample}.coverview
	#/home/programs/CoverView-1.4.3/coverview -i /home/diagnostics/pipelines/Validation/Final_Output/$i/$i".sorted.bam" -b /home/diagnostics/pipelines/Validation/bedfiles/ALLPanel_Bed_Coordinates_hg38_140524_format_sortd.bed -c /home/programs/CoverView-1.4.3/config/config.txt -o /home/diagnostics/pipelines/Validation/Final_Output/$i/${i}.coverview
	#python3 "/home/diagnostics/pipelines/Validation/scripts/coverview.py" /home/diagnostics/pipelines/Validation/Final_Output/$i/${i}.coverview_regions.txt /home/diagnostics/pipelines/Validation/Final_Output/$i/${i}.coverview_regions.csv

done
