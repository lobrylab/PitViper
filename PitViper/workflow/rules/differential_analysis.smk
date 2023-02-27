



rule DESeq2_counts:
    input:
        counts=rules.counts_filtering.output.raw_filtered_counts#config['count_table_file']
    output:
        deseq2_out="results/{token}/DESeq2/{treatment}_vs_{control}/{treatment}_vs_{control}_DESeq2_table.txt"
    params:
        treatment="{treatment}",
        control="{control}",
        design=config['tsv_file']
    log:
        "logs/{token}/DESeq2/{treatment}_vs_{control}.log"
    benchmark:
            "benchmarks/{token}/DESeq2/{treatment}_vs_{control}.tsv"
    script:
        "../../workflow/scripts/DESeq2_analysis.R"



rule directional_scoring_method_method:
    """Implementation of directional scoring method for screened elements prioritization based on sgRNAs filtering."""
    input:
        count_table=rules.counts_filtering.output.normalized_filtered_counts,#config['normalized_count_table'],
        deseq2_table=rules.DESeq2_counts.output.deseq2_out
    output:
        down_elements="results/{token}/directional_scoring_method/{treatment}_vs_{control}/{treatment}_vs_{control}_down-elements_directional_scoring_method.txt",
        up_elements="results/{token}/directional_scoring_method/{treatment}_vs_{control}/{treatment}_vs_{control}_up-elements_directional_scoring_method.txt",
        all_elements="results/{token}/directional_scoring_method/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_directional_scoring_method.txt"
    params:
        treatment="{treatment}",
        control="{control}",
        param_1 = config["directional_scoring_method_fdr_threshold"],
        param_2 = config["directional_scoring_method_log2_threshold"],
        param_3 = config["directional_scoring_method_guides_threshold"]
    log:
        "logs/{token}/directional_scoring_method/{treatment}_vs_{control}.log"
    benchmark:
        "benchmarks/{token}/directional_scoring_method/{treatment}_vs_{control}.tsv"
    script:
        "../../workflow/scripts/directional_scoring_method.R"
