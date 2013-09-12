# Jake Yeung
# 9 Sept 2013
# plot_inclusion_exclusion.R

library(ggplot2)

mydat <- read.table('G:/jyeung/projects/alternative_splicing/output/miso_outputs/mark_rubin/comparison_all_samples/tophat_7520_C3_PE_recurrence.distcalcs.interpolated.txt', header=TRUE, sep='\t')
# mydat <- read.table('G:/jyeung/projects/alternative_splicing/output/miso_outputs/mark_rubin/comparison_all_samples/tophat_8740_PE_recurrence.distcalcs.interpolated.txt', header=TRUE, sep='\t')
str(mydat)

# mydat.inclusion <- mydat[which(mydat$inclusion_or_exclusion == 'inclusion'), ]
qplot(psi, inclusion_or_exclusion, data=mydat, color=group) + scale_shape(solid=FALSE) + geom_jitter(position = position_jitter(height = .03)) + geom_vline(aes(xintercept=0.5), linetype="dashed", size=0.5) + xlab('normalized_psi')