log <- file(snakemake@log[[1]], open="wt")
sink(log, append=TRUE)
sink(log, append=TRUE, type="message")

library(dplyr)
library(readr)
library(fgsea)
library(DESeq2)
library(tibble)
library(matrixStats)


# Handle arguments.
cts_file <- snakemake@input[[1]]
res_file <- snakemake@input[[2]]
treatment <- snakemake@params[1]
baseline <- snakemake@params[2]


# Process counts file (must be tab delimited and have a 'sgRNA' column).
cts <- read.csv(cts_file, sep="\t", row.names="sgRNA", check.names=FALSE)
cts <- cbind(cts[ , grepl( c(baseline) , names( cts ) ) ], cts[ , grepl( c(treatment) , names( cts ) ) ])

columns <- colnames(cts)

control_columns <- columns[grepl(c(baseline), columns)]
treatment_columns <- columns[grepl(c(treatment), columns)]

print(control_columns)
print(treatment_columns)

treatment_df <- cts %>% select(treatment_columns)
control_df <- cts %>% select(control_columns)

table <- data.frame(mean_treatment=rowMeans(treatment_df[,]),
           mean_control=rowMeans(control_df[,]),
           sd_treatment = rowSds(as.matrix(treatment_df[,])),
           sd_control = rowSds(as.matrix(control_df[,]))) %>%
  mutate(rank_metric = (mean_treatment - mean_control) / (sd_treatment + sd_control)) %>%
  arrange(rank_metric) %>%
  rownames_to_column("sgRNA") %>%
  filter(!is.na(rank_metric)) %>%
  select(rank_metric, sgRNA)

ranking <- setNames(table$rank_metric, table$sgRNA)

print(head(ranking))

# Open and read DESeq2 output file.
# res <- read_delim(res_file, "\t", escape_double = FALSE, trim_ws = TRUE)


# Rank sgRNAs using p-value and log2Fold-Change.
# res <- res %>% mutate(rank_metric = log2FoldChange) %>% filter(!is.na(rank_metric)) %>% select(rank_metric, sgRNA)
# ranking <- setNames(res$rank_metric, res$sgRNA)



# Create a set composed of all screened elements.
count_table <- read.table(cts_file, header = TRUE)
gene.sgrna <- count_table[c('sgRNA', 'Gene')]
pathways <- unstack(gene.sgrna)

# Use fgsea library for gsea-like element prioritization.
fgseaRes <- fgsea(pathways = pathways,
                  stats    = ranking,
                  minSize  = 1,
                  maxSize  = Inf)


# Define leadingEdge as a string.
fgseaRes <- fgseaRes %>% select(-leadingEdge)


# Retrieve results.
sign.down.res <- fgseaRes[ES < 0 & pval < 0.05] %>% arrange(pval)
sign.up.res <- fgseaRes[ES > 0 & pval < 0.05] %>% arrange(pval)


write.table(data.frame(sign.down.res), snakemake@output[[1]], quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)
write.table(data.frame(sign.up.res), snakemake@output[[2]], quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)
write.table(data.frame(fgseaRes), snakemake@output[[3]], quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)
