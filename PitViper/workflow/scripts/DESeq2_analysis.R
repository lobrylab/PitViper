log <- file(snakemake@log[[1]], open="wt")
sink(log, append=TRUE)
sink(log, append=TRUE, type="message")


library(dplyr)
library(readr)
library(DESeq2)


# Handle the arguments.
cts_file <- snakemake@input[[1]]
treatment <- snakemake@params[1]
baseline <- snakemake@params[2]


# Process counts file (must be tab delimited and have a 'sgRNA' column).
cts <- read.csv(cts_file, sep="\t", row.names="sgRNA")
cts <- cbind(cts[ , grepl( c(baseline) , names( cts ) ) ], cts[ , grepl( c(treatment) , names( cts ) ) ])


# Create metadata dataframe.
baseline_rep <- unlist(replicate(n = 3, expr = baseline))
treatment_rep <- unlist(replicate(n = 3, expr = treatment))
conditions <- c(baseline_rep, treatment_rep)

cts.names <- colnames(cts)
coldata <- data.frame(row.names = cts.names, condition = conditions)
cts <- cts[, rownames(coldata)]
cts <- na.omit(cts)
cts <- cts[, rownames(coldata)]


# Use DESeq2 for differential analysis.
name_res = paste('condition', treatment, "vs", baseline, sep = "_")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

dds$condition <- factor(dds$condition, levels = c(baseline, treatment))
dds <- DESeq(dds)
print(resultsNames(dds))
res <- results(dds, name=name_res)
res <- data.frame(res)
res <- tibble::rownames_to_column(res, "sgRNA")

# Write output
write_tsv(data.frame(res), snakemake@output[[1]], quote = FALSE)
