con <- file(snakemake@log[[1]])
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")


print(snakemake@input[[1]])
print(snakemake@params[[1]])
print(snakemake@params[[2]])
print(snakemake@params[[3]])

### Libraries
library(readr)
library(dplyr)
library(CRISPhieRmix)

set.seed(123)

### Data
res <- read_delim(snakemake@input[[1]],
                          "\t", escape_double = FALSE, trim_ws = TRUE)
cts <- read_delim(snakemake@input[[2]],
                          "\t", escape_double = FALSE, trim_ws = TRUE)

# Annotate DESeq2 data with elements name
genes <- cts %>% select(Gene, sgRNA)
all_count.DESeq2 <- merge(x = genes, y = res, by = "sgRNA")


if (snakemake@params[[3]] != "") {
  # Open negatives controls guides file
  neg_controls_file <- snakemake@params[[3]]
  # Retrieve neg controls guides ID
  neg_controls_ids <- read.csv(neg_controls_file, header = F, col.names = "id")$id
  # Retrieve neg controls guides log2FoldChange
  neg_ctrl_guides <- which(all_count.DESeq2$sgRNA %in% neg_controls_ids)
  neg_ctrl_guides_log2fc <- all_count.DESeq2 %>% filter(sgRNA %in% neg_controls_ids) %>% pull(log2FoldChange)

  # Retrieve non controls guides log2FoldChanges
  all_count.DESeq2 <- all_count.DESeq2 %>% filter(!sgRNA %in% neg_controls_ids) %>% filter(!is.na(log2FoldChange))
  geneIds = all_count.DESeq2$Gene
  geneIds = factor(geneIds, levels = unique(geneIds))
} else {
  # Retrieve non controls guides log2FoldChanges
  all_count.DESeq2 <- all_count.DESeq2 %>% filter(!is.na(log2FoldChange))
  geneIds = all_count.DESeq2$Gene
  geneIds = factor(geneIds, levels = unique(geneIds))
}


### Compute mean log2FoldChange of Top-3 sgRNAs per target.
meanlog2FoldChanges <- all_count.DESeq2 %>%
  group_by(Gene) %>%
  slice_max(order_by = abs(log2FoldChange), n = 3) %>%
  summarise(top_3_mean_log2FoldChange = mean(log2FoldChange))


### Run CRISPhieRmix
if (snakemake@params[[3]] != "") {
  log2fcCRISPhieRmixFit = CRISPhieRmix::CRISPhieRmix(x = all_count.DESeq2$log2FoldChange, geneIds = geneIds,  negCtrl = neg_ctrl_guides_log2fc)
} else {
  log2fcCRISPhieRmixFit = CRISPhieRmix::CRISPhieRmix(x = all_count.DESeq2$log2FoldChange, geneIds = geneIds)
}
print(log2fcCRISPhieRmixFit)
log2fcCRISPhieRmixScores = data.frame(gene = log2fcCRISPhieRmixFit$genes, locfdr = log2fcCRISPhieRmixFit$locfdr, FDR = log2fcCRISPhieRmixFit$FDR)
log2fcCRISPhieRmixScores$locfdr[which(log2fcCRISPhieRmixScores$FDR < 0)] = 0
log2fcCRISPhieRmixScores = log2fcCRISPhieRmixScores[order(log2fcCRISPhieRmixScores$locfdr, decreasing = FALSE), ]

crisphiermix_res <- merge(x = log2fcCRISPhieRmixScores, y = meanlog2FoldChanges, by.x = "gene", by.y = "Gene")

write.csv(crisphiermix_res, file = snakemake@output[[1]], row.names = FALSE, quote = FALSE)
