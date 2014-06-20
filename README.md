alternative-splicing
====================

This repository contains a number of useful scripts to analyze alternative splicing. Most of these scripts were used to process outputs of MISO. Therefore, these scripts are most useful if you have outputs from MISO to process.

# Brief highlights of this repository:
* Identify differentially spliced events between two groups of samples using t-test.
* Create bed files and fasta files from a list of MISO events.
* Generate heatmaps of PSI values either in pair-wise comparison or group comparison.
* Wrapper function to perform ANCHOR on exons to identify exons that may mediate protein-protein interactions.
* Wrapper function to perform MEME and TOMTOM analysis near genic regions of cassette exons to identify enriched motifs near alternatively spliced events and associate RNA binding proteins to these motifs.
* Perform conservation analysis of enriched motifs by calculating GERP scores of sites contributing to each motif.
* Visualize distribution of sites contributing to each enriched motif.
* Visualize overlap of alternatively spliced events in a Venn diagram.

# Installation

## Packages and libraries required

### Python
* Python 2.6 or greater
* SciPy
* NumPy
### R
* R 3.0.0 or greater
* ggplot2
* extrafont
* grid
* RColorBrewer
* gplots

## Cloning repository
cd dir_to_place_repository
git clone https://github.com/jakeyeung/alternative-splicing.git

# Example workflows

