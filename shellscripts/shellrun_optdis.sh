#!/bin/bash

cd /Data/jyeung/projects/alternative_splicing/git/alternative-splicing/R

copyfiles(){
		cp /Data/jyeung/projects/alternative_splicing/output/tables/$1 /Data/jyeung/projects/alternative_splicing/output/optdis_outputs/$2/activity_gene/MarkerDiscovery/takeda/subtype/Training.txt
		cp /Data/jyeung/projects/alternative_splicing/output/tables/$1 /Data/jyeung/projects/alternative_splicing/output/optdis_outputs/$2/activity_gene/MarkerDiscovery/takeda/subtype/Test.txt
}

fname=gene_expression_data_forOptDis_samplenames_standardized.txt
outputfolder=gene_exprs_only
echo 'Creating folders...'
timeout 30 Rscript run_optdis.R $fname $outputfolder
echo 'Copying files to new folders...'
copyfiles $fname $outputfolder
echo 'Running OptDis'
Rscript run_optdis.R $fname $outputfolder
echo 'Done for gene exprs only'

fname2=si_and_gene_expression_data_forOptDis_withsamplenames_standardized.txt
outputfolder2=si_and_gene_exprs
echo 'Creating folders...'
timeout 30 Rscript run_optdis.R $fname2 $outputfolder2
echo 'Copying files to new folders...'
copyfiles $fname2 $outputfolder2
echo 'Running OptDis'
Rscript run_optdis.R $fname2 $outputfolder2
echo 'Done for si and gene exprs'

fname3=si_and_gene_expression_data_probes_forOptDis_withsamplenames_standardized.txt
outputfolder3=probe_and_gene_exprs
echo 'Creating folders...'
timeout 30 Rscript run_optdis.R si_and_gene_expression_data_probes_forOptDis_withsamplenames_standardized.txt probe_and_gene_exprs
echo 'Copying files to new folders...'
copyfiles $fname3 $outputfolder3
echo 'Running OptDis'
Rscript run_optdis.R $fname3 $outputfolder3
echo 'Done for probe and gene exprs'