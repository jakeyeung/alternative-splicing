#!/bin/bash
#May 29 2014
#meme_run_pipeline_run_psp.sh
#Second step of pipeline: run psp

unshuffled_dir=fasta
shuffled_dir=fasta_shuffled
maindir=/Data/jyeung/projects/alternative_splicing/output/constitutive_exons
#length, mynumb and output_dirname should match meme_run_pipeline_create_background_fastas.sh
mynumb=$1
dirname_suffix=$2
fasta_dir_prefix=$3
psp_dir=$4
fasta_dir=$fasta_dir_prefix/$unshuffled_dir
length=100
constitutives_fasta_dir=$maindir/fasta_files_"$length"bp_"$mynumb"_"$dirname_suffix"

main_dir=/Data/jyeung/projects/alternative_splicing/output/constitutive_exons
# set options
minw=5
maxw=9
# Get fasta from bed
genome_fasta_file=/Data/hg18/Chrs/hg18.fa
# Create fasta files
# pos parameters
mkdir $constitutives_fasta_dir
mkdir $constitutives_fasta_dir/$unshuffled_dir
mkdir $constitutives_fasta_dir/$shuffled_dir

# set dirs
#fasta_dir=/Data/jyeung/projects/alternative_splicing/output/miso_outputs/xenographs_331_331R_hg19_v2/SE.hg19.gff3/filtered_batch/filtered_psi_025/fasta_files_100bp/$unshuffled_dir
bg_dir=$constitutives_fasta_dir/$unshuffled_dir
#psp_dir=/Data/jyeung/projects/alternative_splicing/output/miso_outputs/xenographs_331_331R_hg19_v2/SE.hg19.gff3/filtered_batch/filtered_psi_025/motif_outputs/psp_output/psp_files_100bp_xenograft_rerun_"$mynumb"

# Run PSP before doing meme.

mkdir $psp_dir

for f in `ls $fasta_dir`
do
	echo "Generating PSP for: $f"
	psp-gen -pos $fasta_dir/$f -neg $bg_dir/$f -minw $minw -maxw $maxw -alpha DNA > $psp_dir/$f&
done
wait
echo "All PSPs generated. Files found in: $psp_dir"

