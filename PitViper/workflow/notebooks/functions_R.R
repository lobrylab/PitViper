library(dplyr)
library(RobustRankAggreg)
library(yaml)
library(readr)
library(ggplot2)
library(tidyr)
library(venn)


RobustRankAggregate <- function(ranks) {
    if ("mle_rank" %in% colnames(ranks)) {
        mle <- ranks %>%
            select(id, mle_rank) %>%
            arrange(mle_rank) %>% pull(id)
    }
    
    if ("rra_rank" %in% colnames(ranks)) {
        rra <- ranks %>%
            select(id, rra_rank) %>%
            arrange(rra_rank) %>% pull(id)
    }

    if ("bagel_rank" %in% colnames(ranks)) {
    bagel <- ranks %>%
        select(id, bagel_rank) %>%
        arrange(bagel_rank) %>% pull(id)
    }

    if ("in_house_rank" %in% colnames(ranks)) {
        in_house <- ranks %>%
            select(id, in_house_rank) %>%
            arrange(in_house_rank) %>% pull(id)
    }

    if ("gsea_rank" %in% colnames(ranks)) {
        gsea <- ranks %>%
            select(id, gsea_rank) %>%
            arrange(gsea_rank) %>% pull(id)
    }

    glist <- list(mle, rra, bagel, in_house, gsea)
    glist

    res <- aggregateRanks(glist = glist)

    res <- res %>% mutate(Rank=rank(Score, ties.method= "min"))

    rownames(res) <- NULL

    return(res)
}


mean_counts_by_condition <- function(conditions, gene) {
    conditions <- unlist(strsplit(conditions, ","))

    order.dict <- tsv_file %>%
        filter(condition %in% conditions) %>%
        distinct(condition, replicate)


    replicates <- tsv_file %>%
        filter(condition %in% conditions) %>%
        pull(replicate)


    cts.vertival <- cts_file %>%
        filter(Gene == gene) %>%
        select(Gene, sgRNA, replicates) %>%
        gather(key = "replicate", value = "count", -Gene, -sgRNA)

    merged <- merge(x = cts.vertival, y = order.dict, by = "replicate")

    merged$condition = factor(merged$condition, levels = conditions)


    plot <- merged %>%
        group_by(sgRNA, condition) %>%
        summarize(mean = mean(count)) %>%
        ggplot(aes(x = condition, y = mean, group = sgRNA)) +
            geom_line(aes(color=sgRNA)) +
            geom_point(aes(color=sgRNA)) +
            theme_classic() +
            ggtitle(paste0(gene, ": average number of reads"))
    return(plot)
}


venn_diagram <- function(occ_df, treatment, control) {
    venn(occ_df, ilabels = FALSE, zcolor = "style", ilcs= 1, sncs = 1, borders = FALSE, box = FALSE)
}


get_data <- function(token) {
    config <- paste0("./config/",token,".yaml")
    content = read_yaml(config)
    tsv_file = read_delim(content$tsv_file, "\t", escape_double = FALSE, trim_ws = TRUE)
    cts_file = read_delim(content$normalized_count_table, "\t", escape_double = FALSE, trim_ws = TRUE)

    return(list("tsv" = tsv_file, "cts" = cts_file))
}
