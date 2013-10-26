# Jake Yeung
# 9 Sept 2013
# plot_inclusion_exclusion.R


# Libs --------------------------------------------------------------------


library(ggplot2)



# Functions ---------------------------------------------------------------

plot_dots <- function(df, label, outfile){
    myplot <- 
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
    return(myplot)
    # ggsave(file=outfile, width=2, height=2, scale=3)
}

plot_dots2 <- function(df){
    myplot <- 
        ggplot(df, aes(x=psi, y=iso_status)) +
        geom_point(size=5, position=position_jitter(height=0.08), alpha=0.3) +
        geom_vline(aes(xintercept=0.5), linetype='dashed', size=0.5) +
        geom_vline(aes(xintercept=0), size=0.5) +
        geom_vline(aes(xintercept=1), size=0.5) +
        labs(x=expression(paste('PSI' ['scaled'])),
             y='') +
        theme_bw()
    return(myplot)
}

add_isostatus_remove_pca_nepc <- function(df, sample_label, incl_or_excl){
    # Remove rows containing group PC and NEPC.
    # Rename group sample to sample_label.
    df <- subset(df, group == 'sample' & inclusion_or_exclusion == incl_or_excl)
    df$group <- rep(sample_label, nrow(df))
    iso_status <- apply(df, 1, function(x){
        if(x[3] == 'inclusion'){
            Status <- paste(x[2])
            
        } else if(x[3] == 'exclusion'){
            Status <- paste(x[2])
        } else{
            stop('Neither inclusion or exclusion.')
        }
    })
    df <- cbind(df, iso_status)
    return(df)
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

# Adjust data frames by removing PCa and NEPC mean, adding an iso-status
incl_or_excl <- 'exclusion'
nepc_sample <- add_isostatus_remove_pca_nepc(nepc_sample, 
                                             'NEPC Sample',
                                             incl_or_excl)
pca_sample <- add_isostatus_remove_pca_nepc(pca_sample, 
                                            'PCa Sample',
                                            incl_or_excl)
nepc_hybrid_sample <- add_isostatus_remove_pca_nepc(nepc_hybrid_sample, 
                                                    'NEPC-PCa Mixed',
                                                    incl_or_excl)

# Combine into one dataframe.
master_df <- do.call(rbind, list(nepc_sample, pca_sample, nepc_hybrid_sample))
master_df$iso_status

# Reorder factors
master_df$iso_status <- factor(master_df$iso_status, 
                    levels(factor(master_df$iso_status))[c(3, 2, 1)])
master_df$iso_status
# Plot
gg <- plot_dots2(master_df)

# gg_nepc_sample <- plot_dots(nepc_sample, 'NEPC Sample', nepc_outpath)
# gg_pca_sample <- plot_dots(pca_sample, 'PCa Sample', pca_outpath)
# gg_nepc_hybrid_sample <- plot_dots(nepc_hybrid_sample, 'NEPC-PCa Hybrid Sample', nepc_hybrid_outpath)
