alternative-splicing
====================

This repository contains a number of useful scripts to analyze alternative splicing. Most of these scripts were used to process outputs of MISO. Therefore, these scripts are most useful if you have outputs from MISO to process.

# Brief highlights of this repository, with example workflows:
* [Identify differentially spliced events between two groups of samples using t-test.](#ttest)
* [Create bed files and fasta files from a list of MISO events.](#bedfasta)
* [Generate heatmaps of PSI values either in pair-wise comparison or group comparison.](#heatmaps)
* [Translate nucleotide sequences of cassette exons (or constitutive exons) to amino acid sequences](#createproteinfiles)
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
* [Bio](http://biopython.org/wiki/Main_Page) (for UniProt/SwissProt annotations)

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

### Procedure:
1. [Perform t-test using `t_test_miso_output.py`.](#t_test_miso_output)
2. [Adjust p-values for multiple-test correction (Benjamini-Hochberg method) using `adjust_pvalues.R`](#adjust_pvals)
3. [Append gene names to output file](#append_genenames)

<a name="t_test_miso_output"/>
## 1. Perform t-test on MISO outputs

Inputs for `t_test_miso_output.py`

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

Positional arguments for `append_genenames.py`:

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

### Workflow
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

### Procedure
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

<a name="indexexons"/>
## Index ENSEMBL annotations

Script: `database_scripts/index_exon_info_to_pkl.py`

Takes as input `.gtf` file containing gene annotations from ENSEMBL and stores data as a python dictionary (`.pkl` file)

Positional arguments:

1. `.gtf` file of gene annotations from [ENSEMBL](ftp://ftp.ensembl.org/pub/).
2. Output for pickled python dictionary.

<a name="indexuniprot"/>
## Index UniProt annotations

Script: `database_scripts/index_unitprot_db.py`

Positional arguments:

1. Output pickle path

Other required arguments:
* `-f`: Path for `.dat` file containing uniprot annotations from [UniProt website](http://www.uniprot.org/downloads).

<a name="createproteinfiles"/>
## Translate nucleotide sequences to amino acid sequences
Translate nucleotide sequences to amino acid sequences using `database_scripts/create_dna_protein_summary_file.py`.
Creates a 'protein summary' file that contains the amino acid sequences used for running ANCHOR. The 'protein summary' file also contains information regarding gene names, reading frames and exon numbers for quality control and sanity check purposes.

Option flags:
* `-c --constitutive_exons: use this if your nucleotide sequences come from constitutive exons rather than cassette exons. --help for more information.

Positional arguments:

1. Ensembl dictionary in pkl format (created from [alternative_splicing_scripts/database_scripts/index_exon_info_to_pkl.py](#indexexons)]
2. Fasta files (nucleotides)
3. 1 | 2 | 3. 2 extracts cassette exon. 1 extracts upstream exon. 3 extracts downstream exon.
4. Output path

<a name="runanchor"/>
## Run ANCHOR on amino acid sequences.

Option flags for `motif_scripts/run_anchor_batch.py`:
* `-h, --help for more information`

Positional arguments:

1. Protein summary file from `database_scripts/create_dna_protein_summary_file.py`
2. Output directory
3. Summary output file, placed inside output directory.

Example workflow for running ANCHOR on both inclusion and exclusion exons, see: `example_workflows/run_anchor.sh`

<a name="plotanchor"/>
## Compare ANCHOR results between two sets of exons.

Script: `alternative_splicing_scripts/motif_scripts/plot_anchor_results_comparisons.py`

Option flags:
* `-h, --help for details`

Usage:
`plot_anchor_results_comparisons.py anchor_results1.txt anchor_results2.txt`

<a name="uniprot"/>
## Get UniProt annotations from amino acid sequences

### Procedure
1. Index `.gtf` file and `.dat` file from [ENSEMBL](#indexexons) and [UniProt](#indexuniprot), respectively.
2. [Create protein summary file](#createproteinfiles).
3. [Annotate protein summary file with UniProt annotations](#annotateuniprot).

<a name="annotateuniprot"/>

Script: `annotate_genes_with_swissprot.py`

Outputs: creates a file of amino acids for each exon (if it found a match to the `.gtf` file) and adds UniProt annotations. This allows one to understand the potential functions of protein segments encoded in exons.

Positional arguments:

1. Protein summary file, created from [`database_scripts/create_dna_protein_summary_file.py`](#createproteinfiles).
2. Output file of proteins with UniProt/SwissProt annotations. 

Required flags:
* -f Indexed file of UniProt database (i.e. output from [`database_scripts/index_unitprot_db.py`](#indexuniprot)).

<a name="motifs"/>
## Run motif discovery and identify enriched motifs associated with motifs of RNA binding proteins

Bash script that glues a number of python scripts together: `example_workflows/meme_run_pipeline.full_pipeline.sh`

Outputs: MEME results from the 14 inputs (see Positional Argument 3) and TOMTOM results summarizing the 14 inputs comparing discovered motifs with experimentally derived motifs from [RNA binding proteins](http://cisbp-rna.ccbr.utoronto.ca/). The results can be found in `main_dir/motif_outputs/meme_run1/summary/tomtom.summary` by default.

Positional arguments:

1. Random seed - set to 0 (I found the choice of random seed no longer matters, since we do not randomly take a subset of rows, rather take the entire set of constitutive exons).
2. Evalue - motif with a minimum E-value before MEME stops searching for additional motifs. Play around with this value (10^-5, 10^-10). **Note**: I found that just because MEME stops searching for additional motifs because the motif it just found was above minimum E-value, it does not mean that there does not exist any more motifs wtih E-values below the threshold. Therefore, sometimes setting E-value to something lower (e.g. 10^1) and then manually filtering out motifs that match the E-values may be less prone to false negatives.
3. Main directory containing 3 directories: `fasta_files_100bp/fasta`, `fasta_files_100bp/fasta_shuffled` and `motif_outputs`.
	* `fasta_files_100bp/fasta` contains 14 `.fasta` files (i.e. inputs to MEME discovery):
		* exon_1_inclusion.fasta, exon_1_exclusion.fasta
		* exon_2_inclusion.fasta, exon_2_exclusion.fasta
		* exon_3_inclusion.fasta, exon_3_exclusion.fasta
		* intron_1_5p_inclusion.fasta, intron_1_5p_inclusion.fasta
		* intron_1_3p_inclusion.fasta, intron_1_3p_inclusion.fasta
		* intron_2_5p_inclusion.fasta, intron_2_5p_inclusion.fasta
		* intron_3_3p_inclusion.fasta, intron_3_3p_inclusion.fasta
	* `fasta_files_100bp/fasta_shuffled` contains the same \*.fasta files as in `fasta_files_100bp/fasta` but shuffled using MEME tool `fasta-shuffle-letters`. Names are identical.
	* `motif_outputs`: does not need anything inside. Output from MEME will be found here in directory named `meme_run_seed_$MYSEED_evalue_$MYEVALUE`.
4. Path to MEME database. Format should conform with [MEME](http://meme.nbcr.net/meme/doc/meme-format.html).

<a name="gerp"/>
## Get evolutionary conservation of discovered motifs from MEME

Bash script gluing python scripts together: `example_workflows/get_gerp_scores.sh`

Outputs: Creates a MEME summary file that incorporates GERP conservation scores with discovered motifs.

Positional arguments:

1. Main directory, full directory of motif outputs (e.g. if my MEME outputs are in `main_dir/motif_outputs/meme_run1`) then main directory is `main_dir/motif_outputs`)
2. Name of directory for motif outputs (e.g., if MEME outputs are in `main_dir/motif_outputs/meme_run1`, then set to `meme_run1`)
3. Filename containing TOMTOM summary (e.g., by default it would be the file `tomtom.summary` found in `main_dir/motif_outputs/meme_run1/summary`
4. Directory containing GERP base-wise conservation scores (i.e. \*.maf.rates for chr1, chr2, chr3... chr22, chrX, chrY). They can be downloaded [here](http://mendel.stanford.edu/SidowLab/downloads/gerp/)(**6.3GB!**)

