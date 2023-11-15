rule ssrea_ranking:
    """Ranking for SSREA method."""
    input:
        count_table=rules.counts_filtering.output.normalized_filtered_counts,
        deseq2_table=rules.DESeq2_counts.output.deseq2_out,
        design=config["tsv_file"]
    output:
        "results/{token}/SSREA/{treatment}_vs_{control}/{treatment}_vs_{control}_ranking_SSREA.txt"
    params:
        treatment="{treatment}",
        control="{control}",
        method=config["ssrea_ranking_method"]
    log:
        "logs/{token}/SSREA/ranking_{treatment}_vs_{control}.log"
    message:
        "Preprocessing for SSREA method. From {wildcards.treatment} vs {wildcards.control}."
    script:
        "../../workflow/scripts/ssrea_ranking.py"


rule ssrea:
    """SSREA method for screened elements prioritization."""
    input:
        ranking=rules.ssrea_ranking.output,
        cts=rules.counts_filtering.output.normalized_filtered_counts,
    output:
        "results/{token}/SSREA/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_SSREA.txt"
    log:
        "logs/{token}/SSREA/{treatment}_vs_{control}.log"
    benchmark:
        "benchmarks/{token}/SSREA/{treatment}_vs_{control}.tsv"
    message:
        "SSREA method for screened elements prioritization. From {wildcards.treatment} vs {wildcards.control}."
    script:
        "../../workflow/scripts/ssrea.R"
