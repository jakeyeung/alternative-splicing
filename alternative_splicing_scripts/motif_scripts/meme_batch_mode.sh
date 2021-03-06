#!/bin/bash
# Argument = -m oops|zoops|anr -n <nmotifs> -w <minw> -W <maxw> -e <evalthreshold (for MEME)> -q <qvalue_threshold (for ToMTom)	> -d <database> fasta_dir output_dir

usage()
{
cat << EOF
usage: $0 options

This script runs MEME and TomTom in batch mode.

OPTIONS:
   -h      Show this message
   -m    Distribution of motifs (oops|zoops|anr), default zoops
   -n Number of motifs to find, default 10.
   -w   minimum width of a motif, default 7.
   -W   max width of a motif, default 8.
   -e   E-value threshold for MEME motif discovery, default none.
   -q   Q-value threshold for TomTom motif comparison, default 1.
   -d   Database to match discoverd motifs, default is ray2013_rbp database.
   -M   Maximum size, default 200000
   -p   psp file, from psp-gen output. Corrects against background.
EOF
}

batch_meme()
{
# Run meme on all fasta files in directory.

# Try to make directory for output
mkdir $output_dir
echo $fasta_dir
cd $fasta_dir
# Intialize count
count=0
for f in `ls *.fasta`
do
	echo "Fasta file: $f, PSP file: $psp_dir/$f"
	# Get filename without .fasta extension
	basename="${f%%.*}"
	# Only put -evt option if it is not empty.
	if [ -z "$evt" ]
	then 
		#meme $f -mod $mod -dna -nmotifs $nmotifs -psp $psp_dir/$f -minw $minw -maxw $maxw -maxsize $maxsize -o $output_dir/$basename &
		meme $f -mod $mod -dna -nmotifs $nmotifs -psp $psp_dir/$f -minw $minw -maxw $maxw -maxsize $maxsize -o $output_dir/$basename &
		let count+=1
		[[ $((count%7)) -eq 0 ]] && wait
		
	else
		#meme $f -mod $mod -dna -nmotifs $nmotifs -psp $psp_dir/$f -minw $minw -maxw $maxw -maxsize $maxsize -evt $evt -o $output_dir/$basename &
		meme $f -mod $mod -dna -nmotifs $nmotifs -psp $psp_dir/$f -minw $minw -maxw $maxw -maxsize $maxsize -evt $evt -o $output_dir/$basename &
		let count+=1
		[[ $((count%7)) -eq 0 ]] && wait
	fi
done
wait
echo "All jobs looped for MEME"
}

batch_tomtom()
{
cd $fasta_dir
# Run TomTom on all meme files (run after batch_meme)
for f in `ls *.fasta`
do
	basename="${f%%.*}"
	tomtom -o $output_dir/$basename/rbp_matches -thresh $thresh $output_dir/$basename/meme.txt $db &
done
wait
echo "All jobs looped for TOMTOM"
}

# Get optional arguments...
mod=zoops
nmotifs=10
minw=7
maxw=8
evt=
thresh=1
maxsize=200000
db="/home/jyeung/meme/db/motif_databases/ray2013_rbp.meme"
while getopts “hm:n:w:W:q:d:e:M:p:” OPTION
do
     case $OPTION in
         m)
             mod=$OPTARG
             ;;
         n)
             nmotifs=$OPTARG
             ;;
         w)
             minw=$OPTARG
             ;;
         W)
             maxw=$OPTARG
             ;;
		 q)
			 thresh=$OPTARG
			 ;;
		 d)
			 db=$OPTARG
			 ;;
		 e)
			 evt=$OPTARG
			 ;;
		 M)
			 maxsize=$OPTARG
			 ;;
		 p)
			 psp_dir=$OPTARG
			 ;;
         ?)
             usage
             exit
             ;;
     esac
done

# Get positional arguments...
fasta_dir=${@:$OPTIND:1}
output_dir=${@:OPTIND+1:1}

# Tell user what options are used...
echo "Fasta dir: $fasta_dir"
echo "Output dir: $output_dir"
echo "Distribution of motifs: $mod"
echo "Number of motifs to discover: $nmotifs"
echo "Minimum width of motifs: $minw"
echo "Maximum width of motifs: $maxw"
echo "E-value threshold for MEME: $evt"
echo "Q value threshold for TomTom motif comparison: $thresh"
echo "Database used for TomTom motif comparison: $db"
echo "Maximum size: $maxsize"

# Discover motifs with MEME
# Run meme on all fasta files in directory.

batch_meme
batch_tomtom

