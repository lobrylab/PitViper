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

# Parameters
neg_controls_file <- snakemake@params[[3]]
screen_type <- snakemake@params[[4]]
mu <- as.numeric(snakemake@params[[5]])
bimodal <- snakemake@params[[6]]
top <- as.numeric(snakemake@params[[7]])

if (neg_controls_file != "") {
  # Retrieve neg controls guides ID
  neg_controls_ids <- read.csv(neg_controls_file, header = F, col.names = "id")$id

  # Retrieve neg controls guides log2FoldChange
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

if (top >= 1) {
  ### Compute mean log2FoldChange of Top-3 sgRNAs per target.
  meanlog2FoldChanges <- all_count.DESeq2 %>%
    group_by(Gene) %>%
    slice_max(order_by = abs(log2FoldChange), n = top) %>%
    summarise(mean_log2FoldChange = mean(log2FoldChange))
} else {
  meanlog2FoldChanges <- all_count.DESeq2 %>%
    group_by(Gene) %>%
    slice_max(order_by = abs(log2FoldChange), prop = top) %>%
    summarise(mean_log2FoldChange = mean(log2FoldChange))
}


### Run CRISPhieRmix
if (neg_controls_file != "") {
  log2fcCRISPhieRmixFit = CRISPhieRmix::CRISPhieRmix(x = all_count.DESeq2$log2FoldChange, 
                                                     geneIds = geneIds,
                                                     negCtrl = neg_ctrl_guides_log2fc, 
                                                     screenType = screen_type,
                                                     mu = mu, 
                                                     BIMODAL = bimodal,
                                                     max_iter = 1000,
                                                     tol = 1e-12)
} else {
  log2fcCRISPhieRmixFit = CRISPhieRmix::CRISPhieRmix(x = all_count.DESeq2$log2FoldChange,
                                                     geneIds = geneIds,
                                                     screenType = screen_type,
                                                     mu = mu,
                                                     BIMODAL = bimodal,
                                                     max_iter = 1000,
                                                     tol = 1e-12
                                                     )
}


### Create results table
log2fcCRISPhieRmixScores = data.frame(gene = log2fcCRISPhieRmixFit$genes, locfdr = log2fcCRISPhieRmixFit$locfdr, FDR = log2fcCRISPhieRmixFit$FDR)

### Replace negative FDR values by 0
log2fcCRISPhieRmixScores$locfdr[which(log2fcCRISPhieRmixScores$FDR < 0)] = 0

### Merge results with mean log2FoldChanges
crisphiermix_res <- merge(x = log2fcCRISPhieRmixScores, y = meanlog2FoldChanges, by.x = "gene", by.y = "Gene")

### Write results
write.csv(crisphiermix_res, file = snakemake@output[[1]], row.names = FALSE, quote = FALSE)
