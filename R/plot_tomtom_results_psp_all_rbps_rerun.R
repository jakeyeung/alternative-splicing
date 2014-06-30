# Jake Yeung
# June 10 2014
# plot_tomtom_results_psp_all_rbps_rerun.R
# Plots tomtom results after rerunning results with larger delta Psi cutoff


# Libraries ---------------------------------------------------------------

library(ggplot2)
library("extrafont")    # for Arial font
library(grid)

# LoadFile ----------------------------------------------------------------

tomtom_path <- 'G:/jyeung/projects/ALTERN~1/output/MISO_O~1/XENOGR~1/SEHG19~1.GFF/FILTER~1/FILTER~2/MOTIF_~1/MI949A~1.000/fasta/summary/tomtom.filtered.summary.forplot.txt'
#tomtom_path <- 'G:/jyeung/projects/alternative_splicing/output/xenograft_SE_miso/meme_tomtom_outputs/minw_5_maxw_9_all_rnabps_100bp_psp_all_rbps/fasta_unshuffled/summary/tomtom_concatenated_pcbp3_removed.txt'
tomtom_data <- read.table(tomtom_path, header=TRUE)

str(tomtom_data)

gene_hits <- levels(tomtom_data$gene)


# Save Plot info ----------------------------------------------------------

fname='G:/jyeung/projects/alternative_splicing/figs/tomtom_results_xenograft_arial.rerun.eps'

# Get -log of qvalue ------------------------------------------------------

tomtom_data$pval_log <- -log10(tomtom_data$pval)
tomtom_data$text_pval_log <- -log10(tomtom_data$text_pval)

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

#mygenes <- c('PCBP3', 'CELF3', 'NOVA2', 'ELAVL3', 'SRRM4', 'RBM38', 'PTBP2', 'ESRP2')
#tomtom_data  <- tomtom_data[which(tomtom_data$gene %in% mygenes), ]

# Add dummy variable with intron_2_5p and intron_1_3p for visualization purposes
# dummy_row <- data.frame(gene='', region='intron_2_5p', motif_id='intron_2_3p_inclusion:motif_1', n_sites=NA, pval=0, text_pval=0, eval=0, qval=0, incl_or_excl='inclusion', pval_log=0, text_pval_log=0)
# tomtom_data <- rbind(tomtom_data, dummy_row)
# dummy_row <- data.frame(gene='', region='intron_1_3p', motif_id='intron_2_3p_inclusion:motif_1', n_sites=NA, pval=0, text_pval=0, eval=0, qval=0, incl_or_excl='inclusion', pval_log=0, text_pval_log=0)
# tomtom_data <- rbind(tomtom_data, dummy_row)

# Plot ------------------------------------------------------------------
# plot settings
text_offset <- 0.5    # offset between geom_point and geom_text
pointer_thickness <- 0.75
arrow_thickness <- 0.25    # in cm
# gene modeling drawing settings
y_height <- 0.25
exon1_start <- 0.2
exon1_end <- 1.8
exon2_start <- 3.2
exon2_end <- 4.8
exon3_start <- 6.2
exon3_end <- 7.8
linesize <- 2

loadfonts(device = "postscript")    # make things arial
postscript(fname, height=7, width=6.83, 
           paper='special',
           family='Arial',
           onefile=FALSE,
           horizontal=FALSE)

g <- ggplot(tomtom_data, aes(x=as.numeric(region), y=pval_log, label=gene), legend=FALSE) +
    # add symbols
    geom_point(aes(shape=25, fill=motif_id, size=120)) +
    scale_fill_discrete(guide=FALSE) +
    # label symbosl with text
    geom_text(aes(x=as.numeric(region)+text_offset, y=text_pval_log, size=120, color=motif_id), hjust=0) +
    # draw arrows from text to symbols
    geom_segment(aes(xend=as.numeric(region) + text_offset/10,
                     x=as.numeric(region) + text_offset,
                     yend=as.numeric(pval_log),
                     y=as.numeric(text_pval_log)),
                 colour='grey35',
                 arrow=arrow(length=unit(arrow_thickness,'cm')),
                 size=pointer_thickness,
                 show_guide=FALSE) +
#     geom_segment(data=tomtom_data,
#                  aes(xend=as.numeric(tomtom_data$region),
#                      x=as.numeric(tomtom_data$region) + text_offset,
#                      yend=as.numeric(tomtom_data$pval_log),
#                      y=as.numeric(tomtom_data$text_pval_log)),
#                  colour='grey35',
#                  arrow=arrow(length=unit(arrow_thickness,'cm')),
#                  size=pointer_thickness,
#                  show_guide=FALSE) +
    scale_shape_identity() +
    #scale_size(limits=c(55, 133), range=c(5,10)) +
    scale_y_continuous(limits = c(-0.25, 5.5), 
                       name=expression(paste('-log'[10], ' p-value'))) +
    scale_x_continuous(limits=c(0, 8), breaks=c(1,2,3,4,5,6,7), name='',
                       labels=c("Upstrm\nExon", "Upstrm\nIntron 5'", "Upstrm\nIntron 3'", "Cassette\nExon", "Dwnstrm\nIntron 5'", "Dwnstrm\nIntron 3'", "Dwnstrm\nExon")) + 
#     scale_x_discrete(name='pre-mRNA region',
#                      labels=c("Upstrm Intron 5'", "Upstrm Intron 3'", "Downstrm Intron 5'", "Downstrm Intron 3'")) +
    #labs(size='Number of sites', color='Motif ID') +
    theme_bw() %+replace%
    theme(axis.title=element_text(),
          axis.text.x=element_text(size=10),
          legend.justification=c(1,0), legend.position=c(1,0),
          legend.title=element_text(),
          legend.text=element_text())+
    # put legend at bottom, horizontal
    theme(legend.direction = "horizontal", legend.position = "bottom") + 
    # draw 5' exon
    geom_rect(mapping=aes(xmin=exon1_start, xmax=exon1_end, ymin=-y_height, ymax=y_height), 
              #alpha=0.5, 
              fill='light blue', colour='black', inherit.aes=FALSE, show_guide=FALSE) + 
    # draw cassette exon
    geom_rect(mapping=aes(xmin=exon2_start, xmax=exon2_end, ymin=-y_height, ymax=y_height),
              #alpha=0.5, 
              fill='yellow', colour='black', inherit.aes=FALSE, show_guide=FALSE) +
    # draw 3' exon'
    geom_rect(mapping=aes(xmin=exon3_start, xmax=exon3_end, ymin=-y_height, ymax=y_height),
              #alpha=0.5, 
              fill='light blue', colour='black', inherit.aes=FALSE, show_guide=FALSE) + 
    # draw first intron
    geom_segment(mapping=aes(x=exon1_end, y=0, xend=exon2_start, yend=0), colour='black', size=linesize, show_guide=FALSE) +
    # draw second intron
    geom_segment(mapping=aes(x=exon2_end, y=0, xend=exon3_start, yend=0), colour='black', size=linesize, show_guide=FALSE)
g
dev.off()
print(paste(c('Figure output:', fname)))
