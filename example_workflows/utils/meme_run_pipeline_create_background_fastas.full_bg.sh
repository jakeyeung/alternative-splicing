#!/bin/bash
#May 30 2014
#Run psp-meme pipeline
#instead of taking random subset of rows, take ALL rows (hardcoded)

# Meme library for fasta shuffle letters, in /home/jake
export PATH=${PATH}:${HOME}/meme/bin

#---------
length=100
#name fasta directories
unshuffled_dir=fasta
shuffled_dir=fasta_shuffled
#set random seed
mynumb=$1
#set number of entries in fasta files to extract
#nrows_incl=$2
#nrows_excl=$3
#set suffix filename 
filesuffix=constitutives_pipeline.full_bg.bed
#set directory name for constittuive fasta and beds
#eg: xenograft_pipeline
dirname_suffix=$2
#maindir: /Data/jyeung/projects/alternative_splicing/output/constitutive_exons
#maindir is where you want to output your background bedfiles and fastafiles 
maindir=/Data/jyeung/projects/alternative_splicing/output/constitutive_exons
constitutives_fasta_dir=$maindir/fasta_files_"$length"bp_"$mynumb"_"$dirname_suffix"
output_dir_constitutives_bed=$maindir/bed_files_"$length"bp_"$mynumb"_"$dirname_suffix"
bg_dir=$constitutives_fasta_dir/$unshuffled_dir
dbscriptdir=/Data/jyeung/projects/alternative_splicing/git/alternative-splicing/alternative_splicing_scripts/database_scripts
get_flanking_script=get_flanking_introns_from_constitutive.py
input_dir_constitutives=/Data/jyeung/projects/alternative_splicing/input/buljan_constitutive_exons
input_file=constitutiveExons.txt
#output_dir_constitutives_bed=/Data/jyeung/projects/alternative_splicing/output/constitutive_exons/bed_files_100bp_xenograft_rerun_all
misoscriptdir=/Data/jyeung/projects/alternative_splicing/git/alternative-splicing/alternative_splicing_scripts/miso_scripts
random_script=get_random_subset.py
prefix=constitutiveExons

# Get fasta from bed
genome_fasta_file=/Data/hg18/Chrs/hg18.fa

# set options
minw=5
maxw=9

mkdir $output_dir

for region in upstream downstream exon
do
	python $dbscriptdir/$get_flanking_script --intron_length=$length --include_chr_prefix False --location $region $input_dir_constitutives/$input_file $input_dir_constitutives/"${input_file%.*}".$region."$length"bp."$filesuffix"&
done
wait

mkdir $output_dir_constitutives_bed
#nrows_incl=862
#nrows_excl=504

for region in upstream downstream exon
do
	input_file_bed=$input_dir_constitutives/$prefix.$region."$length"bp."$filesuffix"
	for inclexcl in inclusion exclusion
	do
		if [ "$region" == "exon" ]; then
			# run 3 times, exon1, exon2, exon3
			nrows=56987
			for n in 1 2 3
			do
				echo "Input: $input_dir_constitutives/$prefix.$region.bed"
				echo "Output: $output_dir_constitutives_bed/exon_"$n"_"$inclexcl".bed"
				python $misoscriptdir/$random_script -s $mynumb -n $nrows $input_file_bed $output_dir_constitutives_bed/exon_"$n"_"$inclexcl".bed&
			done
		elif [ "$region" == "upstream" ]; then
			# run 2 times, intron_1_3p and intron_2_3p
			nrows=134545
			for n in 1 2
			do
				echo "Input: $input_dir_constitutives/$prefix.$region.bed"
				echo "Output: $output_dir_constitutives_bed/intron_"$n"_3p_"$inclexcl".bed"
				python $misoscriptdir/$random_script -s $mynumb -n $nrows $input_file_bed $output_dir_constitutives_bed/intron_"$n"_3p_"$inclexcl".bed&
			done
		elif [ "$region" == "downstream" ]; then
			# run 2 times, intron_1_5p and intron_2_5p
			nrows=134545
			for n in 1 2
			do
				echo "Input: $input_dir_constitutives/$prefix.$region.bed"
				echo "Output: $output_dir_constitutives_bed/intron_"$n"_5p_"$inclexcl".bed"
				python $misoscriptdir/$random_script -s $mynumb -n $nrows $input_file_bed $output_dir_constitutives_bed/intron_"$n"_5p_"$inclexcl".bed&
			done
		fi
	done
done
wait

#Create background fasta files

# Create fasta files
# pos parameters
mkdir $constitutives_fasta_dir
mkdir $constitutives_fasta_dir/$unshuffled_dir
mkdir $constitutives_fasta_dir/$shuffled_dir

# set dirs

cd $output_dir_constitutives_bed
for f in `ls *sion.bed`
do
	echo $f
	fastaFromBed -name -s -fi $genome_fasta_file -bed $f -fo $constitutives_fasta_dir/$unshuffled_dir/"${f%.*}.fasta" &
done
wait

