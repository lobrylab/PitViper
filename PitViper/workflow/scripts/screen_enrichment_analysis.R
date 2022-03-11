log <- file(snakemake@log[[1]], open="wt")
sink(log, append=TRUE)
sink(log, append=TRUE, type="message")

library(dplyr)
library(readr)
library(fgsea)
library(DESeq2)
library(tibble)
library(matrixStats)

set.seed(123)

# Handle arguments.
cts_file <- snakemake@input[[1]]
res_file <- snakemake@input[[2]]
treatment <- snakemake@params[1]
baseline <- snakemake@params[2]
method <- snakemake@params[3]
design_file <- snakemake@input[[3]]

# Process counts file (must be tab delimited and have a 'sgRNA' column).
cts <- read.csv(cts_file, sep="\t", row.names="sgRNA", check.names=FALSE)

# Read design file for replicate/condition associations.
design <- read.csv(design_file, sep="\t",  check.names=FALSE)
treatment.name <- design %>% filter(condition == treatment) %>% pull(replicate) %>% unique()
baseline.name <- design %>% filter(condition == baseline) %>% pull(replicate) %>% unique()
cts <- cbind(select(cts, baseline.name), select(cts, treatment.name))

columns <- colnames(cts)

control_columns <- columns[grepl(c(baseline), columns)]
treatment_columns <- columns[grepl(c(treatment), columns)]

n.treatment <- ncol(select(cts, all_of(treatment.name)))
n.baseline <- ncol(select(cts, all_of(baseline.name)))

treatment_df <- cts %>% select(treatment_columns)
control_df <- cts %>% select(control_columns)

if (method == "signal_to_noise" && n.treatment > 2 && n.baseline > 2) {
  print("Signal to Noise")
  table <- data.frame(mean_treatment=rowMeans(treatment_df[,]),
             mean_control=rowMeans(control_df[,]),
             sd_treatment = rowSds(as.matrix(treatment_df[,])),
             sd_control = rowSds(as.matrix(control_df[,]))) %>%
    mutate(rank_metric = (mean_treatment - mean_control) / (sd_treatment + sd_control)) %>%
    arrange(rank_metric) %>%
    rownames_to_column("sgRNA") %>%
    filter(!is.na(rank_metric)) %>%
    select(rank_metric, sgRNA)
} else {
  print("Log2 of classes")
  res <- read_delim(res_file, "\t", escape_double = FALSE, trim_ws = TRUE)
  table <- res %>%
    mutate(rank_metric = log2FoldChange) %>%
    arrange(rank_metric) %>%
    filter(!is.na(rank_metric)) %>%
    select(rank_metric, sgRNA)
}

ranking <- setNames(table$rank_metric, table$sgRNA)

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

print(snakemake@output[[1]])

write.table(data.frame(fgseaRes), snakemake@output[[1]], quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)

