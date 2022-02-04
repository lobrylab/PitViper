con <- file(snakemake@log[[1]])
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")


print(snakemake@input[[1]])
print(snakemake@params[[1]])
print(snakemake@params[[2]])

### Libraries

library(readr)
library(dplyr)
library(CRISPhieRmix)


### Data

res <- read_delim(snakemake@input[[1]],
                          "\t", escape_double = FALSE, trim_ws = TRUE)

cts <- read_delim(snakemake@input[[2]],
                          "\t", escape_double = FALSE, trim_ws = TRUE)


genes <- cts %>% select(Gene, sgRNA)
all_count.DESeq2 <- merge(x = genes, y = res, by = "sgRNA")


log2fc = all_count.DESeq2$log2FoldChange
log2fc.cleaned <- all_count.DESeq2[which(!is.na(all_count.DESeq2$log2FoldChange)),]$log2FoldChange
geneIds = all_count.DESeq2[which(!is.na(all_count.DESeq2$log2FoldChange)),]$Gene
geneIds = factor(geneIds, levels = unique(geneIds))


### Compute mean log2FoldChange of Top-3 sgRNAs per target.

meanlog2FoldChanges <- all_count.DESeq2 %>%
  group_by(Gene) %>%
  slice_max(order_by = abs(log2FoldChange), n = 3) %>%
  summarise(top_3_mean_log2FoldChange = mean(log2FoldChange))



### Run CRISPhieRmix
log2fcCRISPhieRmixFit = CRISPhieRmix::CRISPhieRmix(log2fc.cleaned, geneIds = geneIds)
log2fcCRISPhieRmixScores = data.frame(gene = log2fcCRISPhieRmixFit$genes, locfdr = log2fcCRISPhieRmixFit$locfdr, score = log2fcCRISPhieRmixFit$geneScores, FDR = log2fcCRISPhieRmixFit$FDR)
log2fcCRISPhieRmixScores$locfdr[which(log2fcCRISPhieRmixScores$FDR < 0)] = 0
log2fcCRISPhieRmixScores = log2fcCRISPhieRmixScores[order(log2fcCRISPhieRmixScores$locfdr, decreasing = FALSE), ]

crisphiermix_res <- merge(x = log2fcCRISPhieRmixScores, y = meanlog2FoldChanges, by.x = "gene", by.y = "Gene")

write.csv(crisphiermix_res, file = snakemake@output[[1]], row.names = FALSE, quote = FALSE)
