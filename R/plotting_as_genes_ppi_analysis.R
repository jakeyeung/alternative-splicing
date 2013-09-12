# Jake Yeung
# 9 Sep 2013
# plotting_as_genes_ppi_analysis.R


set.seed(1)
library(ggplot2)
library(plyr)


# Function Defs -----------------------------------------------------------

reduce_rows <- function(main_genes_df, random_genes_df){
    print(colnames(random_genes_df))
    # hope it's > 0
    rows_to_remove <- (nrow(random_genes_df) - nrow(main_genes_df))
    if (rows_to_remove > 0){
        random_rows_indices <- sample(1:rows_to_remove, replace=FALSE)
        random_genes_df.cut <- random_genes_df[-random_rows_indices, , drop=FALSE]
        # colnames(random_genes_df.cut) <- colnames(random_genes_df)
        
    }
    return(random_genes_df.cut)
}

add_rep_col <- function(mydf, new_col_name, new_col_val){
    # Define constants
    dataset <- rep(new_col_val, nrow(mydf))
    mydf.appended <- cbind(mydf, dataset)
    return(mydf.appended)
}

get_mean <- function(jdata, jaes, jfill){
    mean_df <- ddply(jdata, jfill, 
                     my.mean=function(x){mean(x, na.rm=TRUE)}, jaes)
    return(mean_df)
}

plot_histogram <- function(jdata, jaes, jfill){
    m <- ggplot(jdata, 
                 aes(x=jaes, 
                 fill=jfill))
    m + geom_histogram(binwidth=0.2, alpha=.5, position="identity") +
        geom_vline(data=mean_df, aes(xintercept=jaes, colour=jfill),
                   linetype="dashed", size=1)
}

# Read Tables -------------------------------------------------------------

as_results <- read.table('G:/jyeung/projects/alternative_splicing/output/miso_outputs/mark_rubin/comparison_all_samples/ppi_analysis/neighbor_fc_differences_as_genes.txt', 
                         header=TRUE, sep='\t')
random_genes_results <- read.table('G:/jyeung/projects/alternative_splicing/output/miso_outputs/mark_rubin/comparison_all_samples/ppi_analysis/neighbor_fc_differences_random_genes.txt', 
                                   header=TRUE, sep='\t')
str(as_results)
str(random_genes_results)


# Reduce rows from random_gene_results to match as_results size -----------
random_genes_results.cut <- reduce_rows(as_results, random_genes_results)

# Merge as_results and random_genes_results.cut ---------------------------
# so we can ggplot2 easily.
as_results <- add_rep_col(as_results, 'dataset', 'as_gene')
random_genes_results.cut <- add_rep_col(random_genes_results.cut, 
                                        'dataset', 'random_gene')
merged_results <- rbind(as_results, random_genes_results.cut)
ks.test(as_results$number_of_neighbors, random_genes_results.cut$number_of_neighbors)
ks.test(as_results$absolute_gene_exprs_diff_score.degree.normalized., 
       random_genes_results.cut$absolute_gene_exprs_diff_score.degree.normalized.)
ks.test(as_results$as_gene_fold_change, random_genes_results.cut$as_gene_fold_change)
ks.test(as_results$absolute_gene_exprs_diff_score, 
       random_genes_results.cut$absolute_gene_exprs_diff_score)
str(merged_results)

# My Plots: n_neighbors ----------------------------------------------------------------

# Plot my ggplot2
# X logscale
m <- ggplot(merged_results, aes(x=number_of_neighbors, fill=dataset))
m_dens <- m + geom_density(alpha=.3)
# bks <- seq(min(merged_results$number_of_neighbors),max(merged_results$number_of_neighbors))
bks <- c(1, 10, 100, 200, 300)    # My log-scale x scale.
# Find med of each group
cdf <- ddply(merged_results, .(dataset), summarise, n_neighbors.mean=mean(number_of_neighbors))
cdf

# Plot LogScale

m_dens + scale_x_log10(breaks=bks, labels=bks) + geom_vline(data=cdf, 
                                                       aes(xintercept=n_neighbors.mean, 
                                                           colour=dataset),
                                                       linetype="dashed", size=1)

# Plot Linear Scale

m_dens + geom_vline(data=cdf, aes(xintercept=n_neighbors.mean, colour=dataset),linetype="dashed", 
                    size=1)

# Plot Logscale Histogram

m + geom_histogram(binwidth=3, alpha=.5, position="identity") +
    geom_vline(data=cdf, aes(xintercept=n_neighbors.mean, colour=dataset),
               linetype="dashed", size=1)


# My Plots: Gene Exprs ----------------------------------------------------

exprs_plot <- ggplot(merged_results, 
                     aes(x=absolute_gene_exprs_diff_score.degree.normalized., 
                         fill=dataset))
mean_df <- ddply(merged_results, .(dataset), summarise, 
                 abs_gene_exprs_diff.normalized.mean=mean(absolute_gene_exprs_diff_score.degree.normalized., 
                                                          na.rm=TRUE))
exprs_plot + geom_histogram(binwidth=0.2, alpha=.5, position="identity") +
             geom_vline(data=mean_df, aes(xintercept=abs_gene_exprs_diff.normalized.mean, colour=dataset),
                       linetype="dashed", size=1)



# My Plots: Target Gene Exprs ----------------------------------------------------

as_exprs_plot <- ggplot(merged_results, 
                     aes(x=as_gene_fold_change, 
                         fill=dataset))
as_mean_df <- ddply(merged_results, .(dataset), summarise, 
                  as_gene.mean=mean(as_gene_fold_change, na.rm=TRUE))
as_exprs_plot + geom_histogram(binwidth=0.1, alpha=.5, position="identity") +
    geom_vline(data=as_mean_df, aes(xintercept=as_gene.mean, colour=dataset),
               linetype="dashed", size=1)

