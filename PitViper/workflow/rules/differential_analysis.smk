



rule DESeq2_counts:
    input:
        counts=config['count_table_file']
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



rule filtering_method:
    """Implementation of a in-house method for screened elements prioritization based on sgRNAs filtering."""
    input:
        count_table=config['normalized_count_table'],
        deseq2_table=rules.DESeq2_counts.output.deseq2_out
    output:
        down_elements="results/{token}/in_house_method/{treatment}_vs_{control}/{treatment}_vs_{control}_down-elements_in-house.txt",
        up_elements="results/{token}/in_house_method/{treatment}_vs_{control}/{treatment}_vs_{control}_up-elements_in-house.txt",
        all_elements="results/{token}/in_house_method/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_in-house.txt"
    params:
        treatment="{treatment}",
        control="{control}"
    log:
        "logs/{token}/in_house_method/{treatment}_vs_{control}.log"
    benchmark:
        "benchmarks/{token}/in_house_method/{treatment}_vs_{control}.tsv"
    script:
        "../../workflow/scripts/screen_in-house.R"
