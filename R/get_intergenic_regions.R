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
ensembl_or_ucsc <- args[3]
n_per_chr <- 100
ucsc_pos_str_genes <- paste0('/Data/jyeung/projects/alternative_splicing/input',
    '/ucsc_protein_coding_genes/all_pos_strand_genes.txt')
if (ensembl_or_ucsc == 'ensembl' | ensembl_or_ucsc == 'ucsc'){
    print(paste('Found type:', ensembl_or_ucsc))
} else{
    print(paste('Third argument must be either "ensembl" or "ucsc".', 
                ensembl_or_ucsc, 'found.'))
}

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

# Define Chromosomes ------------------------------------------------------

chromosomes <- as.vector(c(seq(1, 22), 'X', 'Y'))
# Add chr prefix to chromosome names
ensembl_or_ucsc='ucsc'
if (ensembl_or_ucsc=='ucsc'){
    i <- 1
    for (mychr in chromosomes){
        chromosomes[i] <- paste0('chr', mychr)
        i <- i + 1
    }
}


# Define Strand -----------------------------------------------------------

strand <- as.vector(c('1'))

# Query biomart or UCSC -----------------------------------------------------------

ucsc_pos_str_genes <- paste0('G:/jyeung/projects/alternative_splicing/input',
                            '/ucsc_protein_coding_genes/all_pos_strand_genes.txt')

# Get all protein coding genes from ENSEMBL
if (ensembl_or_ucsc=='ensembl'){
    strand <- as.vector(c('1'))
    genes <- getBM(attributes=c('chromosome_name', 
                                'start_position', 
                                'end_position',
                                'hgnc_symbol'), 
                   filters=c('chromosome_name', 'strand'),
                   values=list(chromosomes, strand),
                   mart=mart)
} else if (ensembl_or_ucsc=='ucsc'){
    # Get all protein coding genes from UCSC
    genes <- read.table(ucsc_pos_str_genes, header=FALSE, sep='\t', skip=0, 
                        col.names=c('chromosome_name', 
                                    'start_position', 'end_position'))
    # Remove rows that do not have standard chromosome names.
    genes <- subset(genes, chromosome_name %in% chromosomes)
}

# Find complement of protein coding genes ---------------------------------

intergenic_df <- data.frame(matrix(NA, nrow=n_per_chr, ncol=6))
colnames(intergenic_df) <- c('track name=intergenic_regions description="intergenic_regions"',
                             rep('', 5)
count <- 0

# Loop through chromosomes: necessary only for Ensembl, but do it for UCSC for modularity
chromosomes <- c('chr7')
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
        # Add chr prefix only if dealing with ENSEMBL.
        # But either case, add chromosome, start, end, name, score, strand.
        if (ensembl_or_ucsc=='ensembl'){
            chr_start_end_vec <- c(paste0('chr',jchr), start_end_vec, 
                                   'intergene', 0, '+')   
        } else if (ensembl_or_ucsc == 'ucsc'){
            chr_start_end_vec <- c(jchr, start_end_vec, 
                                   'intergene', 0, '+')
        }
        count <- count + 1
        intergenic_df[count, ] <- chr_start_end_vec
    }
}


# Pick random rows from intergenic df -------------------------------------

rand_rows <- sample(1:nrow(intergenic_df), n_rows_to_choose)

intergenic_df.rand_subset <- intergenic_df[rand_rows, ]

# Write to File -----------------------------------------------------------

write.table(intergenic_df.rand_subset, output_fname, quote=FALSE, 
            row.names=FALSE, sep='\t')

print(paste(nrow(intergenic_df.rand_subset), 'rows written to file:', output_fname))