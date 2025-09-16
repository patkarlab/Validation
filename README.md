# Validation

This is a nextflow pipeline for analysing target DNA sequencing data from the LSC WTP samples (CNV myeloid panel)

For running this pipeline, following programs need to be installed and their complete paths need to be added in the params section of the `nextflow.config`.

- adaptors = Fasta file of adapter sequences for trimming
- genome = Genomic fasta file
- samtools = samtools executable path
- bedtools = bedtools executable path
- flt3_bedfile = path to the flt3 bedfile
- ubtf_bedfile =  path to the ubtf bedfile
- site1 = known_polymorphic_sites 1 (Mills_and_1000G_gold_standard.indels)
- site2 = known_polymorphic_sites 2 (dbsnp_138)
- site3 = known_polymorphic_sites 3 (1000G_phase1.snps.high_confidence)
- picard_path = path to the picard.jar file
- GATK38_path = path to the GenomeAnalysisTK-3.8 jar file
- freebayes_path = freebayes executable path 
- platypus_path = path to Platypus.py 
- vardict_path = VarDict executable path
- bcftools_path = bcftools executable path
- strelka_path = path to the strelka bin folder
- lofreq_path = lofreq executable path
- coverview_path = path to the CoverView-1.4.3 folder
- cava_path = path to the CAVA directory
- somaticseq_path = path to the somaticseq_parallel.py
- annovarLatest_path = path to the ANNOVAR folder
- pindel = path to the pindel folder
- get_itd_path = path to the get_itd folder
- java_path = directory containing the java executable
- abra2_path = directory containing the abra jar file
- filt3r_ref = reference fasta for Filt3r
- cnvkit_path = custom CNVkit script path
- vep_script_path = custom VEP script path
- ifcnv = custom ifcnv script path
- deepsomatic = custom deepsomatic script path
- fastp = fastp executable path

## Usage:

1. Keep the `fastq` files into the `sequences/` folder.

2. Change the `samplesheet.csv`. It should have a list of IDs of the samples. 

3. Uncomment the following workflow from `run_nextflow.sh`

```
# For CNV myeloid panel

source activate new_base

nextflow -c  /home/diagnostics/pipelines/Validation/nextflow.config run validation.nf -entry VALIDATION \
--bedfile /home/pipelines/mutation_detector_nextflow/bedfile/CNV_Small_hg19_newmyeloid_sortd \
--cnvkitRef /home/pipelines/MMpanel/scripts/cnvkit_cnvmyeloid/Reference_combpanel.cnn \
--gene_scatter_list /home/pipelines/MMpanel/scripts/cnvkit_cnvmyeloid \
--sequences /home/diagnostics/pipelines/Validation/sequences/ \
--input /home/diagnostics/pipelines/Validation/samplesheet.csv \
-resume -bg

conda deactivate
```

4. `./run_nextflow.sh > script.log`

# NPM1-FLT3 MRD

This workflow is used to detect NPM1 and FLT3 mutations in MRD samples.

For running this workflow, following programs need to be installed and their complete paths need to be added in the params section of the `nextflow.config`.

- get_itd_path = path to the folder containing `getitd.py`
- trimmomatic
- filt3r
- bwa
- samtools = samtools executable path
- bedtools = bedtools executable path
- genome = Genomic fasta file
- bcftools_path = bcftools executable path
- java_path = directory containing the java executable
- varscan_path = path to the VarScan jar file
- annovarLatest_path = path to the ANNOVAR folder
- merge_tsvs_scrip = custom script to merge all the outputs

## Usage:

1. Keep the `fastq` files into the `sequences/` folder.

2. Change the `samplesheet.csv`. It should have a list of IDs of the samples. 

3. Uncomment the following workflow from `run_nextflow.sh`

```
# NPM1-FLT3 amplicon MRD

source activate new_base

nextflow -c /home/diagnostics/pipelines/Validation/nextflow.config run npm1_mrd.nf -entry NPM1 \
--sequences /home/diagnostics/pipelines/Validation/sequences/ \
--input /home/diagnostics/pipelines/Validation/samplesheet.csv \
--bedfile /home/diagnostics/pipelines/Validation/bedfiles/NPM1_FLT3 \
-resume -bg

conda deactivate
```

4. `./run_nextflow.sh > script.log`