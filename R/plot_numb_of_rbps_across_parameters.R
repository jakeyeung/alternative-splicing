# Jake Yeung
# Nov 12 2013
# plot_numb_of_rbps_across_parameters.R
# Reads output from parse_match_rbpid_output.py and
# ggplot2s the number of RBPs across parameter space.


# import_libs -------------------------------------------------------------

library(ggplot2)


# Read table --------------------------------------------------------------

summary_path <- paste0('G:/jyeung/projects/alternative_splicing/output',
    '/motif_outputs/meme_motifs/parameter_space_exploration',
    '/matched_rbpid_summary.txt')
mydf <- read.table(summary_path, header=TRUE, sep='\t')


# Plot --------------------------------------------------------------------

myplot <- ggplot(mydf, 
                 aes(x=MEME_and_filter_parameters, 
                     y=Number_of_RBPs,
                     fill=Shuffled_or_Not)) +
    geom_bar(stat='identity', position='dodge') +
    facet_grid(Location_around_cassette ~ .)
myplot