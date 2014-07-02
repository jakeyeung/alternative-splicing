#!/bin/bash
# June 3 2014
# Get gerp scores from meme motif. Runs two python scripts.
# Define constants
# June 3 2014: Rerun with different thresholds
# Generates both null and interesting.
main_dir=$1
#main_dir="/Data/jyeung/projects/alternative_splicing/output/miso_outputs/rubin_takeda_pooled/SE.hg19.gff3/t_test_results/filter030/motif_outputs"
mydir=$2
#mydir="minw_5_maxw_9_all_rnabps_100bp_psp_all_rbps"
tomtom_filename=$3
#tomtom_filename="tomtom.summary"
gerp_dir=$4
#gerp_dir="/Data/jyeung/projects/alternative_splicing/input/gerp_conservation_score/gerp_scores_by_chr"
scriptdir="/Data/jyeung/projects/alternative_splicing/git/alternative-splicing/alternative_splicing_scripts"
python_script1="motif_scripts/summarize_meme_results.py"
python_script2="motif_scripts/summarize_meme_results_with_gerp_scores.py"
python_appendgene_script="miso_scripts/append_gene_names_to_textfile.py"
annotation_file="/Data/jyeung/projects/alternative_splicing/input/alternative_splicing_annotations/hg19/SE.hg19.gff3"
meme_dir="$main_dir/$mydir/fasta/summary"
meme_files="$main_dir/$mydir/fasta"
tomtom_summary_path="$meme_dir/$tomtom_filename"
#tomtom_summary_path="$meme_dir/tomtom.filtered.summary"
meme_summary_pkl_file1="meme_summary.all.pkl"
meme_summary_out_file1="meme.tomtomfiltered.all.summary"
meme_summary_out_with_gerp="meme.gerp.all.summary"
meme_summary_out_gene_appended="meme.gerp.genename.all.summary"
meme_summary_pkl_file2="gerp_pickle.all.pkl"
gerp_updated_pkl="meme_summary.gerp_updated.all.pkl"
# Run python script 1
echo "Script 1"
echo $meme_dir
echo $tomtom_summary_path
python $scriptdir/$python_script1 -t $tomtom_summary_path -p $meme_summary_pkl_file1 $meme_files $meme_dir/$meme_summary_out_file1
# Run python script 2
echo "Script 2"
python $scriptdir/$python_script2 -p $meme_summary_pkl_file2 -o $gerp_updated_pkl $meme_dir/$meme_summary_pkl_file1 $gerp_dir $meme_dir/$meme_summary_out_with_gerp
# Append gene names to output file
echo "Script 3"
python $scriptdir/$python_appendgene_script -c miso_event $meme_dir/$meme_summary_out_with_gerp $annotation_file $meme_dir/$meme_summary_out_gene_appended
echo "Done"

# Run null analysis on new xenograft.
# March 21 2014
meme_summary_pkl_file1="meme_summary_batchmode.null.pkl"
meme_summary_out_file1="meme_summary_noGERP.null.out"
gerp_dir="/Data/jyeung/projects/alternative_splicing/input/gerp_conservation_score/gerp_scores_by_chr"
meme_summary_out_with_gerp="meme_summary_withGERP.null.out"
meme_summary_out_gene_appended="meme_summary_withGERP.genename_appended.null.out"
gerp_updated_pkl="meme_summary.gerp_updated.null.pkl"
# Run python script 1
python $scriptdir/$python_script1 -t $tomtom_summary_path -n True -A -p $meme_summary_pkl_file1 $meme_files $meme_dir/$meme_summary_out_file1
# Run python script 2
python $scriptdir/$python_script2 -p "gerp_pickle.null.pkl" -o $gerp_updated_pkl $meme_dir/$meme_summary_pkl_file1 $gerp_dir $meme_dir/$meme_summary_out_with_gerp
# Append gene names to output file
python $scriptdir/$python_appendgene_script -c miso_event $meme_dir/$meme_summary_out_with_gerp $annotation_file $meme_dir/$meme_summary_out_gene_appended
