# Plot dist calcs

library(ggplot2)
library(plyr)


# Functions ---------------------------------------------------------------

stderr <- function(x){
    return(sqrt(var(x)/length(x)))
}


# Main --------------------------------------------------------------------


mydf <- read.table('G:/jyeung/projects/alternative_splicing/output/miso_outputs/mark_rubin/comparison_all_samples/hitndrive/distcalcs_summary2.txt', 
                   header=TRUE, sep='\t')
str(mydf)

# Summarize distance means.
# mydfm <- ddply(mydf, .(genelist, phenotype), summarise, avg_dist=mean(dist_mean))
mydfm <- ddply(mydf, .(genelist, phenotype), summarise, avg_dist=mean(dist_mean), 
               sd=sd(dist_mean), n=length(dist_mean))
se <- with(mydfm, sd / sqrt(n))
mydfm <- cbind(mydfm, se)
# Remove BPH because it's not going to be in the paper.
mydfm <- subset(mydfm, phenotype != 'BPH')


# Rename all_aberrant -> "All Differentially Spliced Genes" and drivers -> "Driver Spliced Genes".
levels(mydfm$genelist)[levels(mydfm$genelist) == 'all_aberrant'] <- 'All Spliced Genes'
levels(mydfm$genelist)[levels(mydfm$genelist) == 'drivers'] <- 'Driver Spliced Genes'
# Change NEPC_mixed to NEPC Mixed
levels(mydfm$phenotype)[levels(mydfm$phenotype) == 'NEPC_mixed'] <- 'NEPC-PCa Mixed'

gg_summary <- 
ggplot(mydfm, aes(x=phenotype, y=avg_dist, fill=genelist)) +
    geom_bar(position=position_dodge(), stat='identity') + 
    geom_errorbar(aes(ymin=avg_dist-se, ymax=avg_dist+se),
                  position=position_dodge(0.9), width=0.1,size=0.3) +
    theme_bw() +
    labs(x='Phenotype', y='Distance from PCa', fill='Gene List') +
    opts(legend.direction = "horizontal", legend.position = "bottom")