log <- file(snakemake@log[[1]], open="wt")
sink(log, append=TRUE)
sink(log, append=TRUE, type="message")

library(readr)
library(fgsea)
library(dplyr)

set.seed(123)

table_file <- snakemake@input[[1]]
cts_file <- snakemake@input[[2]]

table <- read.table(table_file, header = TRUE)

ranking <- setNames(table$ranking, table$sgRNA)

# Create a set composed of all screened elements.
count_table <- read.table(cts_file, header = TRUE)
gene.sgrna <- count_table[c('sgRNA', 'Gene')]
pathways <- unstack(gene.sgrna)

# Use fgsea library for gsea-like element prioritization.
fgseaRes <- fgsea(pathways = pathways,
                  stats    = ranking,
                  minSize  = 1,
                  maxSize  = Inf)


# Remove leadingEdge as a string.
fgseaRes <- fgseaRes %>% select(-leadingEdge)

write.table(data.frame(fgseaRes), snakemake@output[[1]], quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)

