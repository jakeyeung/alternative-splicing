# Jake Yeung
# plot_bubble_beltran_vpc_splicing_factors.R
# January 6 2014


# Libraries ---------------------------------------------------------------

library(ggplot2)
library("extrafont")    # for Arial font


# LoadFile ----------------------------------------------------------------

vpc_beltran_path <- "G:/jyeung/projects/alternative_splicing/output/xenograft_exprs/pvals_all_rbps.three_cohorts_combined.txt"
vpc_beltran_data <- read.table(vpc_beltran_path, header=TRUE)



# Plot --------------------------------------------------------------------

str(vpc_beltran_data)
# Filter out REST, it's not an RNA-binding protein
vpc_beltran_data <- vpc_beltran_data[which(vpc_beltran_data$gene != 'REST'), ]


# Filter ------------------------------------------------------------------

fname='G:/jyeung/projects/alternative_splicing/figs/bubble_plot_genes_in_db_all2.pdf'
loadfonts()    # make things arial
pdf(fname, height=11, width=12.83,
   family='Arial', paper='special',
   onefile=FALSE)

ggplot(vpc_beltran_data, aes(x=log10_p_value.beltran, y=log10_p_value.vpc, size=abs_fold_change, color=fc_direction, label=gene), legend=FALSE) +
    geom_point(alpha=0.8) + scale_area(range=c(1,25)) + 
    scale_x_continuous(name=expression(paste('-log'[10], ' BH-Adj P-Value: Beltran Cohort'))) +
    scale_y_continuous(name=expression(paste('-log'[10], ' BH-Adj P-Value: VPC Cohort'))) +
    geom_text(size=vpc_beltran_data$abs_fold_change, color='black', alpha=0.8)+
    labs(size=expression(paste('log'[2], ' FC :331 vs 331R')), color="Fold Change (FC)") + 
    theme_bw(base_size = 25) %+replace%
        theme(legend.justification=c(1,0), legend.position=c(1,0)) 
              #axis.title.x = element_text(size=20), 
              #axis.title.y = element_text(size=20, angle=90),
              #text=element_text(size=20))
dev.off()
print(paste(c('Figure output:', fname)))