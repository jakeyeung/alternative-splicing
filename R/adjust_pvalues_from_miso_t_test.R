# Jake Yeung
# 23 Aug 2013
# Reads miso summary after t_testing and adjusts the p-values.

args <- commandArgs(trailingOnly=TRUE)
filename <- args[1]

# 'G:/jyeung/projects/alternative_splicing/output/miso_outputs/mark_rubin/comparison_all_samples/summary_of_t_tests.txt'
miso_summary <- read.table(filename, 
                           header=TRUE)

pvals <- miso_summary$pval    # Second column.

pvals.adjusted <- p.adjust(pvals, method='BH')

# Create new dataframe with pvals.adjusted as new column
# place new column to right of unadjusted pvals.
output_df <- cbind(event=miso_summary$event, bh_adj_pval=pvals.adjusted, 
                   miso_summary[,-1])

# Save df to a new filename, add .bh_adj.txt to filename.
filename.new <- strsplit(filename, ".txt")
filename.new <- paste0(filename.new, ".bh_adj.txt")
write.table(output_df, file=filename.new, sep='\t')
print(paste('Adjusted', length(pvals.adjusted), 'p-values.'))
print(paste('BH-adjusted file saved to', filename.new))
