#!/bin/bash
# Run anchor through python.
# 10 Jun 2014
# usage: bash run_anchor.sh protein_dir protein_file output_dir output_file
#
protein_dir=$1
protein_file=$2
output_dir=$3
output_file=$4

export PATH=$PATH://home/jyeung/ANCHOR/
export ANCHOR_PATH=/home/jyeung/ANCHOR/
scriptpath=/Data/jyeung/projects/alternative_splicing/git/alternative-splicing/alternative_splicing_scripts/motif_scripts/run_anchor_batch.py
protein_dir=/Data/jyeung/projects/alternative_splicing/output/constitutive_exons/protein_fasta_files_100bp_smaller
output_dir=/Data/jyeung/projects/alternative_splicing/output/constitutive_exons/protein_fasta_files_100bp_smaller/anchor_results
for inclexcl in inclusion exclusion
do
	python $scriptpath $protein_dir/exon_2_"$inclexcl".protein.summary $output_dir anchor.$inclexcl.summary
done
echo "DONE"
