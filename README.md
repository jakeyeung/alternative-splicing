alternative-splicing
====================

This repository contains a number of useful scripts to analyze alternative splicing. Most of these scripts were used to process outputs of MISO. Therefore, these scripts are most useful if you have outputs from MISO to process.

# Brief highlights of this repository, with example workflows:
* [Identify differentially spliced events between two groups of samples using t-test.](#ttest)
* [Create bed files and fasta files from a list of MISO events.](#bedfasta)
* [Generate heatmaps of PSI values either in pair-wise comparison or group comparison.](#heatmaps)
* [Wrapper function to perform ANCHOR on exons to identify exons that may mediate protein-protein interactions.](#anchor)
* [Wrapper function to perform MEME and TOMTOM analysis near genic regions of cassette exons to identify enriched motifs near alternatively spliced events and associate RNA binding proteins to these motifs.](#motifs)
* [Perform conservation analysis of enriched motifs by calculating GERP scores of sites contributing to each motif.](#gerp)
* [Visualize distribution of sites contributing to each enriched motif.](#motifdistribution)
* [Visualize overlap of alternatively spliced events in a Venn diagram.](#venndiagram)

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
```
cd dir_to_place_repository
git clone https://github.com/jakeyeung/alternative-splicing.git
```

<a name="ttest"/>
## T-test to find differentially spliced events between two groups of samples

## Typical procedure:
1. [Perform t-test using `t_test_miso_output.py`.](#t_test_miso_output)
2. [Adjust p-values for multiple-test correction (Benjamini-Hochberg method) using `adjust_pvalues.R`](adjust_pvals)
3. [Append gene names to output file](append_genenames)

<a name="t_test_miso_output">
## 1. Perform t-test on MISO outputs

### Inputs for `t_test_miso_output.py`
* `--group1_file`: textfile containing sample names in group 1
* `--group2_file`: textfile containing sample names in group 2
* `--main_dir`: directory containing sample names (should match group 1 and group 2). Inside each directory contains MISO raw outputs (MISO events catalogued by chromosomes) 
* `--output_directory`: directory for output
* `--output_filename`: name of output file
* `--min_counts: minimum junction read counts in a sample to be considered into a t-test. This value depends on the depth of your RNA-Seq data. Best practices suggest `min_counts` of 10 for a depth in the order of 50 million mapped read pairs

*Example*
```
cd alternative_splicing_scripts/miso_scripts
python t_test_miso_output.py --group1_file /path/to/group_1_samples.txt \ 
							 --group2_file /path/to/group_2_samples.txt \
							 --main_dir /path/containing/miso/outputs \
							 --output_directory output/path/file.txt \ 
							 --min_counts 10
```

<a name="adjust_pvals">
## 2. Adjust p-values for multiple-test correction

### Positional arguments for `adjust_pvals.R`
1. Input file
2. Output file

*Example*
```
cd R
Rscript adjust_pvalues.R output_from_t_test.txt pval_adjusted_output_from_t_test.txt
```

<a name="append_genenames">
## 3. Append gene names to output

### Positional arguments for `append_genenames.py` 
1. Input file
2. Annotations path (gff3 file containing MISO IDs and gene names). One can download the gff3 annotation file from the MISO website
3. Output file

*Example*
```
cd alternative_splicing_scripts/miso_scripts
python append_gene_names_to_textfile.py pval_adjusted_file annotations.gff3 output
```
