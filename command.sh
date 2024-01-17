#!/usr/bin/bash

source activate new_base
nextflow -c  /home/diagnostics/pipelines/Validation/nextflow.config run validation.nf -entry VALIDATION \
--bedfile /home/pipelines/mutation_detector_nextflow/bedfile/CNV_Small_hg19_newmyeloid_sortd \
--trans_bedfile /home/pipelines/MMpanel/bedfiles/MMPanel_translocation_sortd \
--cnvkitRef /home/pipelines/MMpanel/scripts/cnvkit_cnvmyeloid/Reference_combpanel.cnn \
--gene_scatter_list /home/pipelines/MMpanel/scripts/cnvkit_cnvmyeloid \
--sequences /home/diagnostics/pipelines/Validation/sequences/ \
--input /home/diagnostics/pipelines/Validation/samplesheet.csv \
-resume -bg
