# Jake Yeung
# plot_parsed_fimo_outputs.R
# 20 Nov 2013

library(ggplot2)


# get_path ----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
mypath <- args[1]
outpath <- args[2]
min_occur <- as.numeric(args[3])

# Read table --------------------------------------------------------------

mytable <- read.table(mypath, header=TRUE, sep='\t')


# Remove motifs that occur less than 1 times. -----------------------------

mytable <- mytable[which(mytable$n_motif_occurences > min_occur), ]


# Get -log of qvalue ------------------------------------------------------

mytable$pval_log <- -log10(mytable$pval)

# Split inclusion and exclusion -------------------------------------------

# Make inclusions positive and exclusions negative on the n_motif_occurences
# y-axis will split positive and negatives.
adjusted_occurs <- matrix(data=NA, nrow=nrow(mytable), ncol=1)
for (i in 1:nrow(mytable)){
    incl_or_excl <- mytable$inclusion_or_exclusion[i]
    n_motif <- mytable$n_motif_occurences[i]
    if (incl_or_excl == 'exclusion'){
        n_motif <- -1 * as.numeric(n_motif)
}
    adjusted_occurs[i, 1] <- as.numeric(n_motif)
}
mytable$n_motif_occurences <- adjusted_occurs


# Multiply adjusted_occurs by 10 ------------------------------------------


# Reorder factors for intron-exon region ----------------------------------

mytable$exon_intron_region <- factor(mytable$exon_intron_region, 
                                     c('exon_1', 'intron_1_5p', 
                                       'intron_1_3p', 'exon_2', 
                                       'intron_2_5p', 'intron_2_3p', 
                                       'exon_3'))

# Plot with ggplot2 -------------------------------------------------------

myplot <- ggplot(mytable, aes(x=exon_intron_region, y=n_motif_occurences, 
                              color=direction, size=pval_log, 
                              label=rbp_name)) +
        geom_text(alpha=1/2) +
        scale_size(range=c(5,10)) +
        guides(colour = guide_legend(override.aes = list(size=10)))
myplot
ggsave(myplot, file=outpath, width=11.5, height=8)
print(paste('File saved to', outpath))
