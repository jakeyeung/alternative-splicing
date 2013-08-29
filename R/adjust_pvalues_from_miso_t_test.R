# Jake Yeung
# 23 Aug 2013
# Reads miso summary after t_testing and adjusts the p-values.


# ARgs --------------------------------------------------------------------


args <- commandArgs(trailingOnly=TRUE)
filename <- args[1]
bh_pval_cutoff <- as.numeric(args[2])
delta_psi_cutoff <- as.numeric(args[3])


# Functions ---------------------------------------------------------------

filter_events <- function(df, bh_cutoff, delta_psi_cutoff){
    # Filter by bh_adj val.
    df <- subset(output_df, bh_adj_pval < bh_cutoff)
    # Filter by delta_psi_med
    group.vec <- as.vector(df$group)
    psi_med.vec <- as.vector(df$psi_median)
    psi_mean1 <- rep('NA', length(group.vec))
    psi_mean2 <- rep('NA', length(group.vec))
    # Turn to string, then split by comma.
    for(i in 1:length(group.vec)){
        group.list <- strsplit(toString(group.vec[i]), ',')[[1]]
        psi_med.list <- strsplit(toString(psi_med.vec[i]), ',')[[1]]
        g1_ind <- which(group.list == "1")
        g2_ind <- which(group.list == "2")
        # Get psi_med and convert to numeric
        psi1 <- sapply(psi_med.list[g1_ind], as.numeric)
        psi2 <- sapply(psi_med.list[g2_ind], as.numeric)
        # Get mean and store to psi_mean
        psi_mean1[i] <- mean(psi1)
        psi_mean2[i] <- mean(psi2)
    }
    # Append to DF
    df <- cbind(df, psi_mean1)
    df <- cbind(df, psi_mean2)
    delta_psi <- as.numeric(psi_mean2) - as.numeric(psi_mean1)
    abs_delta_psi <- abs(delta_psi)
    df <- cbind(df, delta_psi)
    df <- cbind(df, abs_delta_psi)
    df <- subset(df, abs_delta_psi > delta_psi_cutoff)
    return(df)
}


# Main --------------------------------------------------------------------

# 'G:/jyeung/projects/alternative_splicing/output/miso_outputs/mark_rubin/comparison_all_samples/summary_of_t_tests.txt'
miso_summary <- read.table(filename, 
                           header=TRUE)

pvals <- miso_summary$pval    # Second column.

pvals.adjusted <- p.adjust(pvals, method='BH')

# Create new dataframe with pvals.adjusted as new column
# place new column to right of unadjusted pvals.
output_df <- cbind(event=miso_summary$event, bh_adj_pval=pvals.adjusted, 
                   miso_summary[,-1])

# print(str(output_df))

output_df.filtered <- filter_events(output_df, bh_pval_cutoff, delta_psi_cutoff)


# Filter output for BH and delta psi

# Save df to a new filename, add .bh_adj.txt to filename.
filename.new <- strsplit(filename, ".txt")
filename.new <- paste0(filename.new, ".bh_adj_delta_psi.filtered.txt")
write.table(output_df.filtered, file=filename.new, sep='\t', row.names=FALSE)
print(paste('Adjusted', length(pvals.adjusted), 'p-values.'))
print(paste('BH-adjusted file saved to', filename.new))
