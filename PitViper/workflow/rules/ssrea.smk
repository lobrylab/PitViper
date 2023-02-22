



rule ssrea:
    """SSREA method for screened elements prioritization."""
    input:
        count_table=rules.counts_filtering.output.normalized_filtered_counts,#config['normalized_count_table'],
        deseq2_table=rules.DESeq2_counts.output.deseq2_out,
        design=config["tsv_file"]
    output:
        all_elements="results/{token}/SSREA/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_SSREA.txt"
    params:
        treatment="{treatment}",
        control="{control}",
        method=config["ssrea_ranking_method"]
    log:
        "logs/{token}/SSREA/{treatment}_vs_{control}.log"
    benchmark:
        "benchmarks/{token}/SSREA/{treatment}_vs_{control}.tsv"
    script:
        "../../workflow/scripts/screen_enrichment_analysis.R"
