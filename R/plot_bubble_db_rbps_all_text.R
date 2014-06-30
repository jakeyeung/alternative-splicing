# Jake Yeung
# plot_bubble_db_rps.R
# Plot only RBPs in database.
# March 21 2014

# Libraries ---------------------------------------------------------------

library(ggplot2)
library("extrafont")    # for Arial font


# LoadFile ----------------------------------------------------------------

vpc_beltran_path <- "G:/jyeung/projects/alternative_splicing/output/xenograft_exprs/pvals_db_rbps.three_cohorts_combined.DEBUG.txt"
vpc_beltran_data <- read.table(vpc_beltran_path, header=TRUE)



# Plot --------------------------------------------------------------------

str(vpc_beltran_data)


# Filter ------------------------------------------------------------------

fname='G:/jyeung/projects/alternative_splicing/figs/bubble_plot_genes_in_db_only_final_all_annotated.eps'
loadfonts(quiet=TRUE)    # make things arial
postscript(fname, height=8.83, width=6.83,
#postscript(fname, height=6, width=15,    # for poster
           family='Arial', paper='special',
           horizontal=FALSE,
           onefile=FALSE)

g <- ggplot(vpc_beltran_data, aes(x=log10_p_value.beltran, y=log10_p_value.vpc, size=abs_fold_change, fill=fc_direction), legend=FALSE) +
        geom_point(colour='black', shape=21) + scale_area(range=c(1,10)) + 
        geom_text(aes(label=gene),
                  size=4) +
        scale_x_continuous(name=expression(paste('-log'[10], ' BH-Adj P-Value: Beltran Cohort'))) +
        scale_y_continuous(name=expression(paste('-log'[10], ' BH-Adj P-Value: VPC Cohort'))) +
        #geom_text(size=vpc_beltran_data$abs_fold_change, color='black')+
        labs(size=expression(paste('log'[2], ' FC :331R vs 331')), color="Fold Change (FC)") + 
        theme_bw(base_size = 25) %+replace%
            theme(legend.justification=c(1,0), legend.position=c(1,0)) +
        # put legend at bottom, horizontal
        # opts(legend.position = "right")    # for poster
        opts(legend.direction = "horizontal", legend.position = "bottom")
g
dev.off()
print(paste(c('Figure output:', fname)))