# Plots for HIT'nDRIVE
# Combines two scripts, plot_dist_means.R and plot_inclusion_exclusion.R
# and runs it both. Then multiplots the results.
# plots_for_hitndrive.R

rm(list=ls())

# Source Functions --------------------------------------------------------

# Source multiplot function
source(file.path(getwd(), 'Rutilities', 'multiplot.R'))


# Source Files ------------------------------------------------------------

source(file.path(getwd(), 'plot_dist_means.R'))
source(file.path(getwd(), 'plot_inclusion_exclusion.R'))

# Multiplot four ggplot objects

multiplot(gg, gg_summary, cols=2)