log <- file(snakemake@log[[1]], open="wt")
sink(log, append=TRUE)
sink(log, append=TRUE, type="message")

library(dplyr)
library(stringr)
library(readr)


# Handle arguments.
cts_file <- snakemake@input[[1]]
res_file <- snakemake@input[[2]]
treatment <- snakemake@params[1]
baseline <- snakemake@params[2]


# Open and read input files.
cor <- read.csv(cts_file, , sep="\t") %>% select(sgRNA, Gene)
res <- read.csv(res_file, , sep="\t")


# Merge DESeq2 results table with element-to-sgRNA annotation.
table <- merge(x=cor, y=res, by = "sgRNA")


# Elements prioritization.
in.house.res <- table %>%
  mutate(Efficient = ifelse(abs(log2FoldChange) >= 1 & pvalue <= 0.05, log2FoldChange*(1/-log10(pvalue)), 0)) %>%
  filter(!is.na(Efficient)) %>%
  group_by(Gene) %>%
  summarize(up = sum(Efficient>0), down = sum(Efficient<0), n = n(), score = sum(Efficient)*(down/n), prop = (down/n)) %>%
  arrange(score) %>%
  mutate(category = ifelse( down >= 2 & up >= 2, "ambiguous", ifelse(down >= 2 & up < 2 & score < 0, "down", ifelse(up >= 2 & down < 2  & score > 0, "up", "unchanged"))))


# Get results.
in.house.down <- in.house.res %>%
  filter(category == "down")

in.house.up <- in.house.res %>%
  filter(category == "up")


# Write results.
write.table(data.frame(in.house.down), snakemake@output[[1]], quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)
write.table(data.frame(in.house.up), snakemake@output[[2]], quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)
write.table(data.frame(in.house.res), snakemake@output[[3]], quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)
