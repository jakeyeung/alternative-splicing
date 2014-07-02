#!/bin/bash
#meme_run_meme_tomtom.sh
#run meme and tomtom, summarize afterwards

# Runs meme_batch_mode.sh once with proper max and min widths. (5 and 9)
# Queries SRRM4 database with NOVA2, targeted approach not all splice factors.

export PATH=$PATH:/home/jyeung/meme/bin
# Define script/output directories and paths
# Define min and max widths of motifs.
minw=5
maxw=9
# Create directory to store shuffled and unshuffled discovered motifs.
# user changeable: output directory
length=100
#set meme_out_dirsuffix: e.g. xenograft_run
mynumb=$1
meme_out_dirsuffix=$2
#dir containing fasta files for input of meme
#eg: /Data/jyeung/projects/alternative_splicing/output/miso_outputs/xenographs_331_331R_hg19_v2/SE.hg19.gff3/filtered_batch/filtered_psi_025
maindir=$3
#use psp dir from meme_run_pipeline_run_psp.sh
psp_dir=$4
mydir=minw_"$minw"_maxw_"$maxw"_all_rnabps_"$length"bp_psp_"$mynumb"_"$meme_out_dirsuffix"
evalue=$5
# User Changeable: database name...
# db=/Data/jyeung/projects/alternative_splicing/output/motif_outputs/meme_custom_dbs/indrect_direct_rbps_with_SRRM4.meme
db=$6
# without ESRP: db=/home/jyeung/meme/db/motif_databases/top_splice_factors_withSRRM4.meme
# User changeable: qvalue threshold (for matching RBPs)
q_thres=0.15
# Scripts for motif discovery
scriptdir=/Data/jyeung/projects/alternative_splicing/git/alternative-splicing
bashfile=alternative_splicing_scripts/motif_scripts/meme_batch_mode.sh
fastadir=$maindir/fasta_files_"$length"bp
outputdir=$maindir/motif_outputs
# Scripts for comparing shuffled vs original
pythonfile=alternative_splicing_scripts/motif_scripts/parse_shuffled_unshuffled_evalues.py
rfile=R/plot_eval_distribution.R
# Scripts for matching RBPs
non_null_dir=fasta
null_dir=fasta_shuffled
rbp_dir=rbp_matches
pythonfile_filtertomtom=alternative_splicing_scripts/motif_scripts/filter_tomtom_output.py
tomtom_file=tomtom.txt
rbp_db_path=/Data/jyeung/projects/alternative_splicing/input/rna_binding_motifs/Homo_sapiens_2013_09_13_5-44_pm/RBP_Information_all_motifs.txt
output_file=candidate_rbps.txt
strand=True

mkdir $outputdir/$mydir
# Run MEME and TomTom
for null_or_not_dir in $non_null_dir $null_dir
do
	bash $scriptdir/$bashfile -p $psp_dir -w $minw -W $maxw -M 200000 -e $evalue -d $db $fastadir/$null_or_not_dir $outputdir/$mydir/$null_or_not_dir
done
#bash $scriptdir/$bashfile -w $minw -W $maxw -M 200000 -e 0.0000000001 -d $db $fastadir/$null_dir $outputdir/$mydir/$null_dir

# # Run comparison between shuffled and unshuffled.
comparison_dir=$outputdir/$mydir/unshuffled_shuffled_comparison
comparison_file=comparison.txt
mkdir $comparison_dir
python $scriptdir/$pythonfile $outputdir/$mydir/$non_null_dir $outputdir/$mydir/$null_dir $comparison_dir/$comparison_file
# Read comparison file and plot using ggplot2
Rscript $scriptdir/$rfile $comparison_dir/$comparison_file $comparison_dir/comparison_plot.pdf 

# Match RBPID to RBPNAME in batch mode...
for null_or_not_dir in $non_null_dir $null_dir
do
	joutdir=$outputdir/$mydir/$null_or_not_dir
	for f in `ls $joutdir`
	do
		tomtom_dir=$joutdir/$f/$rbp_dir
		tomtom_path=$tomtom_dir/$tomtom_file
		python $scriptdir/$pythonfile_filtertomtom -q $q_thres -s $strand $tomtom_path $tomtom_dir/$output_file
	done
done

# Summarize TomTom results
summary_script=alternative_splicing_scripts/motif_scripts/summarize_tomtom_results.py
for null_or_not_dir in $non_null_dir $null_dir
do
	non_null_output=$outputdir/$mydir/$null_or_not_dir
	mkdir $non_null_output/summary
	python $scriptdir/$summary_script $non_null_output $non_null_output/summary/tomtom.summary
done
