#!/bin/bash
#May 30 2014
#Run psp-meme pipeline for human cohort
#do not take subset of rows, but take ALL. it is hardcoded so incl_rows and excl_rows are no longer needed.

# Meme library for fasta shuffle letters, in /home/jake
export PATH=${PATH}:${HOME}/meme/bin

#define scripts
create_bg_fastas_script=/Data/jyeung/projects/alternative_splicing/git/alternative-splicing/example_workflows/utils/meme_run_pipeline_create_background_fastas.full_bg.sh
run_psp_script=/Data/jyeung/projects/alternative_splicing/git/alternative-splicing/example_workflows/utils/meme_run_pipeline_run_psp.sh
meme_tomtom_script=/Data/jyeung/projects/alternative_splicing/git/alternative-splicing/example_workflows/utils/meme_run_meme_tomtom.user_input_eval.run_shuffled_and_unshuffled.sh

#define inputs for scripts
seed=$1
evalue=$2
dir_suffix=meme_run_seed_"$seed"_evalue_"$evalue"
#main_dir=/Data/jyeung/projects/alternative_splicing/output/miso_outputs/rubin_takeda_pooled/SE.hg19.gff3/t_test_results/filter030
main_dir=$3
fasta_dir=$main_dir/fasta_files_100bp
psp_dir=$main_dir/motif_outputs/psp_output/psp_files_run_"$seed"

#requries two inputs: $1=setseed $2=dirname_suffix
bash $create_bg_fastas_script $seed $dir_suffix
#run psp
bash $run_psp_script $seed $dir_suffix $fasta_dir $psp_dir
#run meme and tomtom
bash $meme_tomtom_script $seed $dir_suffix $main_dir $psp_dir $evalue
