#!/bin/bash
# Argument = -mod oops|zoops|anr -nmotifs <nmotifs> -minw <minw> -maxw <maxw> 

usage()
{
cat << EOF
usage: $0 options

This script runs MEME and TomTom in batch mode.

OPTIONS:
   -h      Show this message
   -mod    Distribution of motifs (oops|zoops|anr)
   -nmotifs Number of motifs to find.
   -minw   minimum width of a motif.
   -maxw   max width of a motif.
EOF
}

mod=
nmotifs=
minw=
maxw=
while getopts “:mod:nmotifs:minw:maxw” OPTION
do
     case $OPTION in
         mod)
             mod=$OPTARG
             ;;
         nmotifs)
             nmotifs=$OPTARG
             ;;
         minw)
             minw=$OPTARG
             ;;
         maxw)
             maxw=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

# if [[ -z $mod ]] || [[ -z $nmotifs ]] || [[ -z $minw ]] || [[ -z $maxw]]
# then
     # usage
     # exit 1
# fi

echo $mod $nmotifs $minw $maxw