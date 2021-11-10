library(dplyr)
library(RobustRankAggreg)

pitviper_R <- function() {
    message("Hello, world!")
}


RobustRankAggregate <- function(df_merged_reduced) {
    mle <- df_merged_reduced %>%
        select(id, mle_rank) %>%
        arrange(mle_rank) %>% pull(id)

    rra <- df_merged_reduced %>%
        select(id, rra_rank) %>%
        arrange(rra_rank) %>% pull(id)

    bagel <- df_merged_reduced %>%
        select(id, bagel_rank) %>%
        arrange(bagel_rank) %>% pull(id)

    in_house <- df_merged_reduced %>%
        select(id, in_house_rank) %>%
        arrange(in_house_rank) %>% pull(id)

    gsea <- df_merged_reduced %>%
        select(id, gsea_rank) %>%
        arrange(gsea_rank) %>% pull(id)

    glist <- list(mle, rra, bagel, in_house, gsea)
    glist

    res <- aggregateRanks(glist = glist)

    res <- res %>% mutate(Rank=rank(Score, ties.method= "min"))

    rownames(res) <- NULL

    return(res)
}
