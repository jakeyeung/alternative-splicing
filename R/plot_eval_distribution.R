# Jake Yeung
# Oct 26 2013
# plot_evalues.R
# From get_e_values_from_meme_output.py -> plot the output in ggplot2 to compare
# evalues from "non-null" vs "null" sequence motifs.


# ImportLibs --------------------------------------------------------------

library(ggplot2)


# Def Consts --------------------------------------------------------------

# myfile <- paste0('G:/jyeung/projects/alternative_splicing/output/',
#                  'motif_outputs/meme_motifs/SE.hg19.gff3/mincount_10/',
#                  'exon_1_exclusion.evalues.parsed')
# myfile <- paste0('G:/jyeung/projects/alternative_splicing/output/',
#                  'motif_outputs/meme_motifs/meme_comparisons/',
#                  'default_comparisons.out')

args <- commandArgs(trailingOnly=TRUE)
myfile <- args[1]
outfile <- args[2]

# ReadFiles ---------------------------------------------------------------

eval_df <- read.table(myfile, header=TRUE, sep='\t')
str(eval_df)

# -log transform ----------------------------------------------------------

transformed_evals <- apply(eval_df, 1, function(x){
    eval <- as.numeric(x[2])
    return(-log10(eval))
    }
)
# Bind to column

eval_df <- cbind(eval_df, transformed_evals)

# Plot --------------------------------------------------------------------

myplot <- ggplot(eval_df, aes(x=transformed_evals)) + 
            geom_density(aes(group=Group, colour=Group, fill=Group), alpha=0.7) + 
            # scale_x_log10() + 
            xlab('-log(E-Value)') + 
            ylab('Density') + 
            theme_bw(30) +
            opts(legend.position="bottom")
ggsave(filename=outfile, plot=myplot, width=11.5, height=8)
