# Jake Yeung
# 9 Sept 2013
# plot_inclusion_exclusion.R


# Libs --------------------------------------------------------------------


library(ggplot2)



# Functions ---------------------------------------------------------------

plot_dots <- function(df, label, outfile){
    ggplot(df, aes(x=psi, y=inclusion_or_exclusion, color=group)) +
        geom_point(size=5, shape=1, position=position_jitter(height=0.08)) +
        geom_vline(aes(xintercept=0.5), linetype='dashed', size=0.5) +
        labs(x=expression(paste('PSI' ['scaled'])), 
             y='Isoform Status', 
             color='Phenotype') +
        scale_color_discrete(name='Phenotype',
                             breaks=c('NEPC', 'PC', 'sample'),
                             labels=c('Avg NEPC', 'Avg PCa', label)) +
        theme_bw() +
        opts(legend.direction = "horizontal", legend.position = "bottom")
    ggsave(file=outfile, width=2, height=2, scale=3)
}


# Set filenames -----------------------------------------------------------


nepc_path <- paste0('G:/jyeung/projects/alternative_splicing/',
                    'output/miso_outputs/mark_rubin/comparison_all_samples/hitndrive/',
                    'tophat_7800_T76_PE.hitndrive.interpolated.txt')
pca_path <- paste0('G:/jyeung/projects/alternative_splicing/',
                   'output/miso_outputs/mark_rubin/comparison_all_samples/hitndrive/',
                   'tophat_97_T_PE.hitndrive.interpolated.txt')
nepc_hybrid_path <- paste0('G:/jyeung/projects/',
                           'alternative_splicing/output/miso_outputs/mark_rubin/',
                           'comparison_all_samples/hitndrive/',
                           'tophat_7520_C3_PE.hitndrive.interpolated.txt')

nepc_outpath <- paste0('G:/jyeung/projects/alternative_splicing/',
    'output/miso_outputs/mark_rubin/comparison_all_samples/hitndrive/',
    'plots/NEPC_sample_dots.pdf')
pca_outpath <- paste0('G:/jyeung/projects/alternative_splicing/',
                   'output/miso_outputs/mark_rubin/comparison_all_samples/hitndrive/',
                   'plots/PCA_sample_dots.pdf')
nepc_hybrid_outpath <- paste0('G:/jyeung/projects/alternative_splicing/',
                              'output/miso_outputs/mark_rubin/comparison_all_samples/hitndrive/',
                              'plots/NEPC_PCa_hybrid_dots.pdf')

nepc_sample <- read.table(nepc_path, header=TRUE, sep='\t')
pca_sample <- read.table(pca_path, header=TRUE, sep='\t')
nepc_hybrid_sample <- read.table(nepc_hybrid_path, header=TRUE, sep='\t')

plot_dots(nepc_sample, 'NEPC Sample', nepc_outpath)
plot_dots(pca_sample, 'PCa Sample', pca_outpath)
plot_dots(nepc_hybrid_sample, 'NEPC-PCa Hybrid Sample', nepc_hybrid_outpath)
