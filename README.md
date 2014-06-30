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

### Computational biology software
* [MEME version 4.9.1 or greater](http://meme.nbcr.net/meme/doc/meme-install.html)
* [ANCHOR](http://anchor.enzim.hu/) 
* [bedtools](http://bedtools.readthedocs.org/en/latest/)

## Cloning repository
```
cd dir_to_place_repository
git clone https://github.com/jakeyeung/alternative-splicing.git
```

<a name="ttest"/>
## T-test to find differentially spliced events between two groups of samples

## Typical procedure:
1. [Perform t-test using `t_test_miso_output.py`.](#t_test_miso_output)
2. [Adjust p-values for multiple-test correction (Benjamini-Hochberg method) using `adjust_pvalues.R`](#adjust_pvals)
3. [Append gene names to output file](#append_genenames)

<a name="t_test_miso_output"/>
## 1. Perform t-test on MISO outputs

### Inputs for `t_test_miso_output.py`
* `--group1_file`: textfile containing sample names in group 1
* `--group2_file`: textfile containing sample names in group 2
* `--main_dir`: directory containing sample names (should match group 1 and group 2). Inside each directory contains MISO raw outputs (MISO events catalogued by chromosomes) 
* `--output_directory`: directory for output
* `--output_filename`: name of output file
* `--min_counts`: minimum junction read counts in a sample to be considered into a t-test. This value depends on the depth of your RNA-Seq data. Best practices suggest `min_counts` of 10 for a depth in the order of 50 million mapped read pairs

*Example*
```
cd alternative_splicing_scripts/miso_scripts
python t_test_miso_output.py --group1_file /path/to/group_1_samples.txt \ 
							 --group2_file /path/to/group_2_samples.txt \
							 --main_dir /path/containing/miso/outputs \
							 --output_directory output/path/file.txt \ 
							 --min_counts 10
```

<a name="adjust_pvals"/>
## 2. Adjust p-values for multiple-test correction

### Positional arguments for `adjust_pvals.R`
1. Input file
2. Output file

*Example*
```
cd R
Rscript adjust_pvalues.R output_from_t_test.txt pval_adjusted_output_from_t_test.txt
```

<a name="append_genenames"/>
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

<a name="bedfasta"/>
## Generate bed and fasta files from MISO outputs
Currently, these scripts only work for MISO outputs for skipped exons/cassette exons.

## Workflow
Scripts used to create bed files:
* `alternative_splicing_scripts/miso_scripts/extract_coordinates_from_miso_summary.py`
* `alternative_splicing_scripts/miso_scripts/split_beds_into_inclusion_exclusion.py`
Scripts and functions used to create fasta files:
* `fastaFromBed`
* `alternative_splicing_scripts/fasta_scripts/remove_chr_from_bed_file.py` (if human genome fasta file does not contain "chr" prefix in front of chromosome locations, you can remove them using this script.)

Quick and dirty bash script example gluing these scripts together: `example_workflows/get_beds_fastas.full_pipeline.sh`

<a name="heatmaps"/>
## Generate heatmaps from MISO outputs.
Works for MISO outputs for any types of alternative splicing. These scripts allow visualization of MISO outputs (e.g. differentially spliced events) in a heatmap. Works for both pairwise comparison as well as groupwise comparison, but reshaping the data into a matrix format differs. Workflows for both types of comparisons are described below.

## Procedure
1. [Reshape MISO output textfile into a matrix suitable for plotting in R.](#reshape)
2. [Read matrix file, plot in R.](#plotheatmap)

<a name="reshape"/>
## Reshaping MISO output depends on whether MISO output was generated from a [pairwise comparison (Bayes factor)](#reshapebf) or [groupwise comparison (t-test)](#reshapettest)
<a name="reshapebf"/>
## 1a. Reshape MISO output to matrix: pairwise comparison (Bayes Factor)
Inputs for alternative_splicing_scripts/miso_scripts/prepare_data_for_clustering_misobf.py:

Arguments from option flags:
* `--sample1_name`: name of sample 1
* `--sample2_name`: name of sample 2

Positional arguments:
1. MISO filtered output for pairwise comparison using Bayes Factor
2. Output file

<a name="plotheatmap"/>
### Plot heatmap using `R/plot_heatmap_psi_values.R`

Example:
`Rscript R/plot_heatmap_psi_values.R reshaped_miso_output.txt myfigure.eps`

* Input file for `R/plot_heatmap_psi_values.R` should be the output after [reshaping](#reshape).

<a name="reshapettest"/>
## 1b. Reshape MISO output to matrix: groupwise comparison (t-test)
Positional arguments:
1. MISO output file generated from [groupwise comparison via t-test](#ttest)
2. Filename containing sample names for group 1
3. Filename containing sample names for group 2
4. Output file

<a name="plotheatmap"/>
Positional arguments for R/plot_heatmap_psi_values.R:
1. MISO output in matrix format [(reshaped)](#reshape)
2. Output EPS file

<a name="anchor"/>
## Run ANCHOR analysis on MISO outputs.
Written for analysis of cassette exons with ANCHOR. ANCHOR predicts protein binding regions within disordered segments.

## Procedure
1. [Translate nucleotide sequences (stored as `.fasta` format) to amino acid sequences.](#createproteinfiles)
2. [Run ANCHOR on amino acid sequences.](#runanchor)
3. [Plot enrichment of cassette exons compared to constitutive exons](#plotanchor) 

<a name="createproteinfiles"/>
## Translate nucleotide sequences to amino acid sequences
Translate nucleotide sequences to amino acid sequences using `alternative_splicing_scripts/database_scripts/create_dna_protein_summary_file.py`.

Option flags:
* `-c --constitutive_exons: use this if your nucleotide sequences come from constitutive exons rather than cassette exons. --help for more information.

Positional arguments:
1. Ensembl dictionary in pkl format (created from [alternative_splicing_scripts/database_scripts/index_exon+info_to_pkl.py](#indexexons)]
2. Fasta files (nucleotides)
3. 1 | 2 | 3. 2 extracts cassette exon. 1 extracts upstream exon. 3 extracts downstream exon.
4. Output path

<a name="runanchor">
## Run ANCHOR on amino acid sequences.

Option flags:
* `-h, --help for more information`

Positional arguments:
1. protein summary file from `create_dna_protein_summary_file.py`
2. Output directory
3. Summary output file, placed inside output directory.

<a name="plotanchor"/>
## Compare ANCHOR results between two sets of exons.

Script: `alternative_splicing_scripts/motif_scripts/plot_anchor_results_comparisons.py`

Option flags:
* `-h, --help for details`

Usage:
`plot_anchor_results_comparisons.py anchor_results1.txt anchor_results2.txt`

