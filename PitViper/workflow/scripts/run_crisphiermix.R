
con <- file(snakemake@log[[1]])
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")


list.of.packages <- c('CRISPhieRmix')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]
if(length(new.packages)) {
  devtools::install_github('timydaley/CRISPhieRmix')
}

print(snakemake@input[[1]])
print(snakemake@params[[1]])
print(snakemake@params[[2]])

### Libraries

library(readr)
library(DESeq2)
library(dplyr)
library(CRISPhieRmix)


### Data

all_count <- read_delim(snakemake@input[[1]],
                        "\t", escape_double = FALSE, trim_ws = TRUE)

n_control <- snakemake@params[[1]]
n_treatment <- snakemake@params[[2]]

control <- replicate(n = n_control, expr = 0)
treatment <- replicate(n = n_treatment, expr = 1)

conditions <- c(control, treatment)


counts = all_count[3:length(all_count)]
colnames(counts)
coldata = data.frame(condition = conditions)
rownames(coldata) = colnames(counts)
all_count.DESeq2 <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

all_count.DESeq2 = DESeq2::DESeq(all_count.DESeq2)
all_count.DESeq2 = DESeq2::results(all_count.DESeq2)

genes <- all_count %>% select(Gene)
all_count.DESeq2 <- cbind(genes, all_count.DESeq2)

log2fc = all_count.DESeq2$log2FoldChange
log2fc.cleaned <- all_count.DESeq2[which(!is.na(all_count.DESeq2$log2FoldChange)),]$log2FoldChange
geneIds = all_count.DESeq2[which(!is.na(all_count.DESeq2$log2FoldChange)),]$Gene
geneIds = factor(geneIds, levels = unique(geneIds))

log2fcCRISPhieRmixFit = CRISPhieRmix::CRISPhieRmix(log2fc.cleaned, geneIds = geneIds, PLOT = TRUE, VERBOSE = TRUE)

print(log2fcCRISPhieRmixFit)

#log2fcCRISPhieRmixScores = data.frame(gene = log2fcCRISPhieRmixFit$genes, FDR = log2fcCRISPhieRmixFit$FDR)
log2fcCRISPhieRmixScores = data.frame(gene = log2fcCRISPhieRmixFit$genes, locfdr = log2fcCRISPhieRmixFit$locfdr, score = log2fcCRISPhieRmixFit$geneScores, FDR = log2fcCRISPhieRmixFit$FDR)
log2fcCRISPhieRmixScores$locfdr[which(log2fcCRISPhieRmixScores$locfdr < 0)] = 0
log2fcCRISPhieRmixScores = log2fcCRISPhieRmixScores[order(log2fcCRISPhieRmixScores$locfdr, decreasing = FALSE), ]

write.csv(log2fcCRISPhieRmixScores, file = snakemake@output[[1]], row.names = FALSE, quote = FALSE)
write.csv(all_count.DESeq2, file = snakemake@output[[2]], row.names = FALSE, quote = FALSE)
