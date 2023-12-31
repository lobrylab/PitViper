log <- file(snakemake@log[[1]], open="wt")
sink(log, append=TRUE)
sink(log, append=TRUE, type="message")


library(dplyr)
library(readr)
library(DESeq2)


# Handle arguments.
cts_file <- snakemake@input[[1]]
treatment <- snakemake@params[1]
baseline <- snakemake@params[2]
design_file <- snakemake@params[[3]]

# Set seed for reproducibility.
set.seed(123)

# Process counts file (must be tab delimited and have a 'sgRNA' column).
cts <- read.csv(cts_file, sep="\t", row.names="sgRNA", check.names=FALSE)
print(colnames(cts))

# Read design file for replicate/condition associations.
design <- read.csv(design_file, sep="\t",  check.names=FALSE)
treatment.name <- design %>% filter(condition == treatment) %>% pull(replicate) %>% unique()
baseline.name <- design %>% filter(condition == baseline) %>% pull(replicate) %>% unique()
cts <- cbind(select(cts, baseline.name), select(cts, treatment.name))
n.treatment <- ncol(select(cts, all_of(treatment.name)))
n.baseline <- ncol(select(cts, all_of(baseline.name)))

# Create metadata dataframe.
baseline_rep <- unlist(replicate(n = n.baseline, expr = baseline))
treatment_rep <- unlist(replicate(n = n.treatment, expr = treatment))
conditions <- c(baseline_rep, treatment_rep)

# Reorder columns to match metadata.
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

# Run DESeq2.
dds$condition <- factor(dds$condition, levels = c(baseline, treatment))
dds <- DESeq(dds)
res <- results(dds, name=name_res)
res <- data.frame(res)
res <- tibble::rownames_to_column(res, "sgRNA")

# Write output
print(snakemake@output[[1]])
write.table(x = data.frame(res), file = snakemake@output[[1]], append = FALSE, sep = "\t", quote = FALSE, row.names = FALSE)
