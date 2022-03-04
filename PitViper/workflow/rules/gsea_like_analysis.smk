



rule gsea_like:
    """Implementation of a GSEA-like method for screened elements prioritization."""
    input:
        count_table=rules.counts_filtering.output.normalized_filtered_counts,#config['normalized_count_table'],
        deseq2_table=rules.DESeq2_counts.output.deseq2_out,
        design=config["tsv_file"]
    output:
        down_elements="results/{token}/GSEA-like/{treatment}_vs_{control}/{treatment}_vs_{control}_down-elements_GSEA-like.txt",
        up_elements="results/{token}/GSEA-like/{treatment}_vs_{control}/{treatment}_vs_{control}_up-elements_GSEA-like.txt",
        all_elements="results/{token}/GSEA-like/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_GSEA-like.txt"
    params:
        treatment="{treatment}",
        control="{control}"
    log:
        "logs/{token}/GSEA-like/{treatment}_vs_{control}.log"
    benchmark:
        "benchmarks/{token}/GSEA-like/{treatment}_vs_{control}.tsv"
    script:
        "../../workflow/scripts/screen_enrichment_analysis.R"
