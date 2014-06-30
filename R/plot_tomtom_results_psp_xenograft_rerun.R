# Jake Yeung
# Jan 9 2014
# plot_tomtom_results.R


# Libraries ---------------------------------------------------------------

library(ggplot2)


# LoadFile ----------------------------------------------------------------

#tomtom_path <- 'G:/jyeung/projects/ALTERN~1/output/MISO_O~1/XENOGR~1/SEHG19~1.GFF/FILTER~1/NUA142~1.BOT/MEME_T~1/XENOGR~4/fasta_unshuffled/summary/tomtom.summary2'
tomtom_path <- 'G:/jyeung/projects/alternative_splicing/output/xenograft_SE_miso/meme_tomtom_outputs/minw_5_maxw_9_all_rnabps_100bp_psp/fasta_unshuffled/tomtom.summary'
tomtom_data <- read.table(tomtom_path, header=TRUE)

str(tomtom_data)

gene_hits <- levels(tomtom_data$gene)

# Get -log of qvalue ------------------------------------------------------

tomtom_data$pval_log <- -log10(tomtom_data$pval)

# Split inclusion and exclusion -------------------------------------------

# Make inclusions positive and exclusions negative on the n_motif_occurences
# y-axis will split positive and negatives.
adjusted_occurs <- matrix(data=NA, nrow=nrow(tomtom_data), ncol=1)
for (i in 1:nrow(tomtom_data)){
    incl_or_excl <- tomtom_data$incl_or_excl[i]
    pval_log <- tomtom_data$pval_log[i]
    if (incl_or_excl == 'exclusion'){
        pval_log <- -1 * as.numeric(pval_log)
    }
    adjusted_occurs[i, 1] <- as.numeric(pval_log)
}
tomtom_data$pval_log <- adjusted_occurs

# Reorder factors for intron-exon region ----------------------------------

tomtom_data$region <- factor(tomtom_data$region, 
                             c('exon_1', 'intron_1_5p', 
                               'intron_1_3p', 'exon_2', 
                               'intron_2_5p', 'intron_2_3p', 
                               'exon_3'))


# Filter genes ------------------------------------------------------------

mygenes <- c('PCBP3', 'CELF3', 'NOVA2', 'ELAVL3', 'SRRM4', 'RBM38', 'PTBP2', 'ESRP2')
tomtom_data  <- tomtom_data[which(tomtom_data$gene %in% mygenes), ]

# Add dummy variable with intron_2_5p and intron_1_3p for visualization purposes
dummy_row <- data.frame(gene='', region='intron_2_5p', motif_id='intron_2_3p_inclusion:motif_1', n_sites=60, pval=0, eval=0, qval=0, incl_or_excl='inclusion', pval_log=0)
tomtom_data <- rbind(tomtom_data, dummy_row)
dummy_row <- data.frame(gene='', region='intron_1_3p', motif_id='intron_2_3p_inclusion:motif_1', n_sites=60, pval=0, eval=0, qval=0, incl_or_excl='inclusion', pval_log=0)
tomtom_data <- rbind(tomtom_data, dummy_row)

# Plot ------------------------------------------------------------------

g <- ggplot(tomtom_data, aes(x=region, y=pval_log, size=n_sites, color=motif_id, label=gene), legend=FALSE) +
    geom_text(alpha=0.8, position=position_jitter(width=0.0, height=0.5))+
    scale_size(range=c(5,10)) +
    scale_y_continuous(limits = c(-4, 4), 
                       name=expression(paste('-log'[10], ' p-value'))) +
    scale_x_discrete(name='pre-mRNA region',
                     labels=c("Upstrm Intron 5'", "Upstrm Intron 3'", "Downstrm Intron 5'", "Downstrm Intron 3'")) +
    labs(size='Number of sites', color='Motif ID') +
    theme_bw() %+replace%
    theme(axis.title=element_text(size=24),
          axis.text=element_text(size=18),
          legend.justification=c(1,0), legend.position=c(1,0))

#axis.title.x = element_text(size=20), 
#axis.title.y = element_text(size=20, angle=90),
#text=element_text(size=20))
y_height <- 0.5
exon1_start <- 0
exon1_end <- 0.9
exon2_start <- 2.1
exon2_end <- 2.9
exon3_start <- 4.1
exon3_end <- 5
linesize <- 2

# Draw 5' exon
g <- g + geom_rect(mapping=aes(xmin=exon1_start, xmax=exon1_end, ymin=-y_height, ymax=y_height), 
                   alpha=0.5, fill='light blue', colour='black', inherit.aes=FALSE, show_guide=FALSE)
# Draw cassette exon
g <- g + geom_rect(mapping=aes(xmin=exon2_start, xmax=exon2_end, ymin=-y_height, ymax=y_height),
                   alpha=0.5, fill='yellow', colour='black', inherit.aes=FALSE, show_guide=FALSE)
# Draw 3' exon
g <- g + geom_rect(mapping=aes(xmin=exon3_start, xmax=exon3_end, ymin=-y_height, ymax=y_height),
                   alpha=0.5, fill='light blue', colour='black', inherit.aes=FALSE, show_guide=FALSE)
# Draw first intron
g <- g + geom_segment(mapping=aes(x=exon1_end, y=0, xend=exon2_start, yend=0), colour='black', size=linesize, show_guide=FALSE)
# Draw second intron
g <- g + geom_segment(mapping=aes(x=exon2_end, y=0, xend=exon3_start, yend=0), colour='black', size=linesize, show_guide=FALSE)

g
#fname='G:/jyeung/projects/alternative_splicing/figs/motif_plot.pdf'
#ggsave(plot=g, file=fname, width=17, height=10, scale=1)