# Jake Yeung
# get_hgnc_symbol_from_coordinates.R
# July 31 2013


# Libraries ---------------------------------------------------------------

library(biomaRt)



# Constants ---------------------------------------------------------------

input_file <- 'G:/jyeung/projects/alternative_splicing/output/miso_outputs/mark_rubin/comparison_4pc_4nepc/pca_four_only_vs_nepc.miso_bf_1_1_10_02_10.filtered'
# input_file <- 'G:/jyeung/projects/alternative_splicing/output/miso_outputs/mark_rubin/comparison_4pc_4nepc/processed_files/pca_four_only_vs_nepc.miso_bf_1_1_10_02_10.bed'
output_file <- 'G:/jyeung/projects/alternative_splicing/output/miso_outputs/mark_rubin/comparison_4pc_4nepc/pca_four_only_vs_nepc.miso_bf_1_1_10_02_10_genenames_appended.filtered'


# LoadTable ---------------------------------------------------------------

miso_file <- read.table(input_file, header=FALSE, sep='\t', 
                        skip=1, stringsAsFactors=FALSE)


# DefineBiomart -----------------------------------------------------------

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")


# GetCoordinates ----------------------------------------------------------

# Replace chr1 with 1.

# chromosome_fixed <- sapply(miso_file[, 1], function(x){
#     x_fixed <- substr(x, 4, nchar(x))
#     return(x_fixed)
# })

regions <- apply(miso_file, 1, function(x){
    # Get regions from id
    miso_id <- x[1]    # $event_name
    # Split by @
    miso_id_parsed <- strsplit(miso_id, "@")[[1]][1]    # Get middle exon.
    # Remove first three chars (chr) and last two chars from string.
    miso_id_rm <- substr(miso_id_parsed, 4, nchar(miso_id_parsed)-2)
    return(miso_id_rm)
})


# Define Attributes -------------------------------------------------------

attributes <- c("wikigene_name","chromosome_name","start_position", 
                "end_position", "strand")
# filters <- c("chromosome_name", "start", "end")
# results <- getBM(attributes=attributes,filters=filters,
#                  values=list(chromosome_fixed[1:10], miso_file[,2][1:10], miso_file[,3][1:10]),
#                  mart=mart)
# str(results)
results <- getBM(attributes=attributes,filters="chromosomal_region",
                 values=regions,mart=mart)
str(results)