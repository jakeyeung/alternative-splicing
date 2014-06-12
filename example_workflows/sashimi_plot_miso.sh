#!/bin/bash
# Sashimi plots
# usage: bash sashimi_plot_miso.sh type_of_splicing coordinate_to_plot ymax
# types of plots: SE.hg19.gff3 MXE.hg19.gff3 etc

type=$1
coords=$2
#ymax=$3
outfile=$3
# Edit settingsfile
settings=/Data/jyeung/projects/alternative_splicing/input/alternative_splicing_annotations/sashimi_plot_settings.txt
sed -i "s/^\(miso_prefix\s*=\s*\/mnt\/enclosure\/jyeung\/miso_outputs\/xenographs_331_331R_hg19_v2\/\s*\).*$/\1$type/" $settings
#sed -i "s/^\(ymax\s*=s*\).*\$/\1$ymax/" $settings

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$ROOTSYS/lib:/mnt/enclosure/mofan/software/ATLAS/ATLAS_LINUX/lib/:/usr/lib64/python2.6/:/usr/local/lib/:/usr/local/lib/libffi-3.0.10/include/:/mnt/enclosure/mofan/software/bcbio-nextgen-0.2/git/bcbb/nextgen/bcbio
annots=/mnt/enclosure/jyeung/miso_annotations/indexed_gff3_hg19/$type
outdir=/mnt/enclosure/jyeung/miso_outputs/xenographs_331_331R_hg19_v2/sashimi_plots

plot.py --plot-event $coords $annots $settings --output-dir $outdir

# Replace colons with underscores, then copy to /Data

data_outdir=/Data/jyeung/projects/alternative_splicing/output/miso_outputs/xenographs_331_331R_hg19_v2/sashimi_plots
cd $outdir
mv "$coords.pdf" $data_outdir/"$outfile"
echo "$coords.pdf moved to $data_outdir/'$outfile'"
