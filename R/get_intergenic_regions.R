# Jake Yeung
# October 15 2013
# Get intergenic regions from biomart.

rm(list=ls())

library(biomaRt)


# Set Seed ----------------------------------------------------------------

set.seed(0)

# Define Biomart Object ---------------------------------------------------

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")


# Set Arguments -----------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
output_fname <- args[1]
n_rows_to_choose <- args[2]
n_per_chr <- 100

# My Functions ------------------------------------------------------------

get_random_intergenic_coordinate <- function(ordered_bm_output, 
                                             flanking_bp=300){
    # Inputs:
    # ordered_bm_output: BM output, ordered by start then end position.
    # flanking_bp: basepair 
    # Outputs:
    # chromosome name, start, end, position.
    # 
    # Intergene start is coding gene end, intergene end is coding gene start.
    # Assumes bm_output is ordered by start position and contains + strands only.
    
    # Set constants
    valid_region <- FALSE    # Initialize
    
    while (valid_region == FALSE){
        # Pick random row
        random_row_i <- sample(2:nrow(ordered_bm_output), 1)
        # Calculate intergene start and end, ASSUMES AN ORDERED DF!
        intergene_start <- ordered_bm_output[(random_row_i-1), 'end_position']
        intergene_end <- ordered_bm_output[(random_row_i), 'start_position']
        # Check if intergene start/end is valid (i.e. end > start)
        if (intergene_end > intergene_start){
            valid_region = TRUE
        }
    }
    # Pick random chromosomal position within intergene region
    # Must subtract flanking_bp to make sure final region does not go outside
    # of intergene region.
    random_region_start <- sample((intergene_start):(intergene_end-flanking_bp), 1)
    # Assumes all
    random_region_end <- random_region_start + flanking_bp
    return(c(random_region_start, random_region_end))
}

# Define Constants --------------------------------------------------------

# n_per_chr <- 100    # Choose 100 regions in each chromosome.
# n_rows_to_choose <- 50    # Pick how many intergenic regions
# output_fname <- paste0("G:/jyeung/projects/alternative_splicing/output/",
#     "miso_outputs/mark_rubin_hg19_v2_rl_insertdist/SE.hg19.gff3/",
#     "t_test_results/bed_files/intergenic_regions.bed")

# Define Chromosomes ------------------------------------------------------

chromosomes <- as.vector(c(seq(1, 22), 'X', 'Y'))
strand <- as.vector(c('1'))

# Query biomart -----------------------------------------------------------
# Get all protein coding genes.
genes <- getBM(attributes=c('chromosome_name', 
                            'start_position', 
                            'end_position'), 
               filters=c('chromosome_name', 'strand'),
               values=list(chromosomes, strand),
               mart=mart)
str(genes)


# Find complement of protein coding genes ---------------------------------

intergenic_df <- data.frame(matrix(NA, nrow=n_per_chr, ncol=6))
colnames(intergenic_df) <- c('track name=intergenic_regions description="intergenic_regions"',
                             '', '', '', '', '')
count <- 0

# Loop through chromosomes
for(jchr in chromosomes){
    df.subset <- subset(genes, chromosome_name == jchr)
    # Order by start, then by end.
    df.subset.sort <- df.subset[with(df.subset, order(start_position, 
                                                      end_position)),]
    # Remove duplicate rows
    df.subset.sort <- unique(df.subset.sort)
    
    # Get intergenic region
    for (i in 1:n_per_chr){
        start_end_vec <- get_random_intergenic_coordinate(df.subset.sort)
        chr_start_end_vec <- c(paste0('chr',jchr), start_end_vec, 
                               'intergene', 0, '+')
        count <- count + 1
        intergenic_df[count, ] <- chr_start_end_vec
    }
}
str(intergenic_df)


# Pick random rows from intergenic df -------------------------------------

rand_rows <- sample(1:nrow(intergenic_df), n_rows_to_choose)

intergenic_df.rand_subset <- intergenic_df[rand_rows, ]

str(intergenic_df.rand_subset)


# Write to File -----------------------------------------------------------

write.table(intergenic_df.rand_subset, output_fname, quote=FALSE, 
            row.names=FALSE, sep='\t')