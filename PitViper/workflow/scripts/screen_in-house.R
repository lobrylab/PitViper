log <- file(snakemake@log[[1]], open="wt")
sink(log, append=TRUE)
sink(log, append=TRUE, type="message")

library(dplyr)
library(stringr)
library(readr)

set.seed(123)

# Handle arguments.
cts_file <- snakemake@input[[1]]
res_file <- snakemake@input[[2]]
treatment <- snakemake@params[1]
baseline <- snakemake@params[2]
inhouse_fdr_threshold <- snakemake@params[3]
inhouse_log2_threshold <- snakemake@params[4]
inhouse_guides_threshold <- snakemake@params[5]


# Open and read input files.
cor <- read.csv(cts_file, , sep="\t") %>% select(sgRNA, Gene)
res <- read.csv(res_file, , sep="\t")


# Merge DESeq2 results table with element-to-sgRNA annotation.
table <- merge(x=cor, y=res, by = "sgRNA")

# Find minimum non-zero p-adjusted value in table.
min.fdr <- table %>% filter(padj > 0) %>% pull(padj) %>% min()

# Replace p-adjusted values equal to 0 by minimal non-zero p-adjusted value.
table <- table %>% mutate(padj = ifelse(padj == 0, min.fdr, padj))

# Elements prioritization.
in.house.res <- table %>%
  mutate(Efficient = ifelse(abs(log2FoldChange) >= inhouse_log2_threshold & padj <= inhouse_fdr_threshold, log2FoldChange*(-log10(padj)), 0)) %>%
  filter(!is.na(Efficient)) %>%
  group_by(Gene) %>%
  summarize(up = sum(Efficient>0), down = sum(Efficient<0), n = n(), scores = list(Efficient)) %>%
  mutate(category = ifelse(down >= inhouse_guides_threshold & up >= inhouse_guides_threshold, "ambiguous",
                    ifelse(down >= inhouse_guides_threshold & up < inhouse_guides_threshold, "down",
                    ifelse(up >= inhouse_guides_threshold & down < inhouse_guides_threshold, "up", "unchanged")))) %>%
  rowwise() %>%
  mutate(score = ifelse(category == "down" | category == "up", sum(scores)/(up + down), 0)) %>%
  arrange(score) %>%
  mutate(scores = 0)

print(in.house.res)

# Get results.
in.house.down <- in.house.res %>%
  filter(category == "down")

in.house.up <- in.house.res %>%
  filter(category == "up")


# Write results.
write.table(data.frame(in.house.down), snakemake@output[[1]], quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)
write.table(data.frame(in.house.up), snakemake@output[[2]], quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)
write.table(data.frame(in.house.res), snakemake@output[[3]], quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)
