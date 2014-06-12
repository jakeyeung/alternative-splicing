#!/bin/bash	
# May 28 2014
# Get bed files and fasta files for 331vs331R: 100bp long.
# User changeable parameters: main_dir and comparison_results

# Meme library for fasta shuffle letters, in /home/jake
export PATH=${PATH}:${HOME}/meme/bin

#main_dir=/Data/jyeung/projects/alternative_splicing/output/miso_outputs/xenographs_331_331R_hg19_v2/SE.hg19.gff3/filtered_batch/filtered_psi_030
main_dir=$1
miso_bf_file=$2
#type=bf OR type=ttest ONLY
type=$3
scriptdir=/Data/jyeung/projects/alternative_splicing/git/alternative-splicing/alternative_splicing_scripts/miso_scripts
extractscript=extract_coordinates_from_miso_summary.py
splitscript=split_beds_into_inclusion_exclusion.py
comparison_results=$main_dir/$miso_bf_file
length=100
bed_dir=$main_dir/bed_files_"$length"bp
mkdir $bed_dir
split_dir=incl_excl_beds
# Create bed files
python $scriptdir/$extractscript -l $length $comparison_results $bed_dir
mkdir $bed_dir/$split_dir
python $scriptdir/$splitscript -t $type $bed_dir $bed_dir/$split_dir $comparison_results

# Create fasta files
# pos parameters
nochr_dir=chr_removed
mkdir $bed_dir/$split_dir/$nochr_dir
fastadir=$main_dir/fasta_files_"$length"bp
unshuffled_dir=fasta
shuffled_dir=fasta_shuffled
mkdir $fastadir
mkdir $fastadir/$unshuffled_dir
mkdir $fastadir/$shuffled_dir
# Get fasta from bed
genome_fasta_file=/mnt/enclosure/mofan/database/HG19/Homo_sapiens.GRCh37.62.dna.chromosome.fa
# Remove chr from bed
remove_chr_script=/Data/jyeung/projects/alternative_splicing/git/alternative-splicing/alternative_splicing_scripts/fasta_scripts/remove_chr_from_bed_file.py
cd $bed_dir/$split_dir
for f in `ls *sion.bed`
do
	python $remove_chr_script $f $bed_dir/$split_dir/$nochr_dir/$f &
done
wait
cd $bed_dir/$split_dir/$nochr_dir
for f in `ls *sion.bed`
do
	fastaFromBed -name -s -fi $genome_fasta_file -bed $f -fo $fastadir/$unshuffled_dir/"${f%.*}.fasta" &
done
wait

# Shuffle fasta files
cd $fastadir/$shuffled_dir
for f in `ls ../$unshuffled_dir/*.fasta`
do
	F=$(basename $f)
	fasta-shuffle-letters < $f > "${F%%.*}".shuffled.fasta
done
