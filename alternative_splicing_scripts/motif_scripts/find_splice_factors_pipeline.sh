#!/bin/bash
# Create custom motif database from a list differentially expressed splice factors 
# (as well as ones Anna recommended).
# Runs FIMO in batchmode. Uses top splice factor motifs as database.
# Parse FIMO outputs and get a summary file from that.
# ggplot2 the summary file.

export PATH=$PATH:/home/jyeung/meme/bin/

usage()
{
cat << EOF
usage: $0 options

This script creates DB of splice factors, runs FIMO from custom DB, then plots results.

OPTIONS:
   -h      Show this message
   -r      Directory containing RBP motif info. Default /Data/jyeung/projects/alternative_splicing/input/rna_binding_motifs/Homo_sapiens_2013_09_13_5-44_pm
   -s      Name of parsed summary file. Default parsed_summary.out
   -t      FDR threshold for FIMO to call a significant "hit". Default 0.15
   -p      Plot file name. Default fimo_summary.pdf
EOF
}

create_custom_db(){
# Create custom motif db
# $1: gene_list=/Data/jyeung/projects/alternative_splicing/output/rubin_gene_exprs_outputs/top_splice_factors_list_plus_indirectdirect_ensemblid.txt
# $2: rbp_motif_dir=/Data/jyeung/projects/alternative_splicing/input/rna_binding_motifs/Homo_sapiens_2013_09_13_5-44_pm
# $3: customdb_output=/Data/jyeung/projects/alternative_splicing/output/motif_outputs/meme_custom_dbs/top_splice_factors_plus_top_indirectdirect.fromscript.meme
buildscript=/Data/jyeung/projects/alternative_splicing/git/alternative-splicing/alternative_splicing_scripts/motif_scripts/build_custom_meme_db.py
python $buildscript $1 $2 $3
}

run_fimo_batch(){
# Run FIMO in batch mode
# $1: outputdir=/Data/jyeung/projects/alternative_splicing/output/motif_outputs/fimo_outputs/mincount_10_fdr_15_DE_splicefactors
# $2: fastadir=/Data/jyeung/projects/alternative_splicing/output/miso_outputs/mark_rubin_hg19_v2_rl_insertdist/SE.hg19.gff3/t_test_results/fasta_files
# $3: thresh=0.15
# $4: customdb_output=/Data/jyeung/projects/alternative_splicing/output/motif_outputs/meme_custom_dbs/top_splice_factors_plus_top_indirectdirect.fromscript.meme
# unshuffleddir=mincount_10
# shuffleddir=mincount_10_shuffledbymeme_commandline
mkdir $1
# non_null_dir=mincount_10
# null_dir=mincount_10_shuffled
for s in $unshuffleddir $shuffleddir
do
	mkdir $1/$s
	for f in `ls $2/$s/*.fasta`
	do
		fasta_path="${f%.*}"
		fasta_basename=$(basename $fasta_path)
		fimo -o $1/$s/$fasta_basename --norc --qv-thresh --thresh $3 $4 $f &
	done
done
wait
}

run_summary_file(){
# Parse FIMO output
# $1: outputdir=/Data/jyeung/projects/alternative_splicing/output/motif_outputs/fimo_outputs/mincount_10_fdr_15_DE_splicefactors
# $2: summaryfile=parsed_summary.out
parsescript=/Data/jyeung/projects/alternative_splicing/git/alternative-splicing/alternative_splicing_scripts/motif_scripts/parse_fimo_outputs.py
python $parsescript $1 $2
}

plot_parsed_output(){
# Plot in R the output
# $1: summaryfile
# $2: outputfile
plotRscript=/Data/jyeung/projects/alternative_splicing/git/alternative-splicing/R/plot_parsed_fimo_outputs.R
Rscript $plotRscript $1 $2
}

# Get optional arguments...
rbp_motif_dir=/Data/jyeung/projects/alternative_splicing/input/rna_binding_motifs/Homo_sapiens_2013_09_13_5-44_pm
summaryfile=parsed_summary.out
unshuffleddir=mincount_10
shuffleddir=mincount_10_shuffledbymeme_commandline
non_null_dir=mincount_10
null_dir=mincount_10_shuffled
thresh=0.15
plotfile=fimo_summary.pdf
while getopts “hr:t:s:p:” OPTION
do
     case $OPTION in
         r)
             rbp_motif_dir=$OPTARG
             ;;
         t)
             thresh=$OPTARG
             ;;
         s)
             summaryfile=$OPTARG
             ;;
         p)
             plotfile=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

# Get positional arguments...
gene_list=${@:$OPTIND:1}
fasta_dir=${@:$OPTIND+1:1}
customdb_output=${@:OPTIND+2:1}
output_dir=${@:OPTIND+3:1}

# Tell user what options/positionals are used...
echo "Optional parameters are set to:"
echo "RBP motif dir: $rbp_motif_dir"
echo "FDR Threshold for FIMO: $thresh"
echo "Positional arguments used:"
echo "Gene list: $gene_list"
echo "Fasta directory: $fasta_dir"
echo "Output directory: $output_dir"
echo "Custom DB Output: $customdb_output"

create_custom_db $gene_list $rbp_motif_dir $customdb_output
run_fimo_batch $output_dir $fasta_dir $thresh $customdb_output
for d in $unshuffleddir $shuffleddir
do
	run_summary_file $output_dir/$d $output_dir/$d/$summaryfile
	plot_parsed_output $output_dir/$d/$summaryfile $output_dir/$d/$plotfile 100
done