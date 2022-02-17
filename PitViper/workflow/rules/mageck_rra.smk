


rule mageck_rra:
    """MAGeCK RRA implementation for gene essentiality analysis.
        input: count table and design matrix
        output: genes and sgrnas summaries
        params: user must choose a normalization method"""
    input:
        count_table=rules.counts_filtering.output.normalized_filtered_counts#config['normalized_count_table']
    output:
        gene_summary = "results/{token}/MAGeCK_RRA/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt",
    params:
        name = "results/{token}/MAGeCK_RRA/{treatment}_vs_{control}/{treatment}_vs_{control}",
        treatment = getTreatmentIds,
        control = getControlIds,
        lfc_method_opt = config['mageck_rra_LFC'],
        ajd_method_opt = config['mageck_rra_adj'],
        rm_zero_threshold_opt = config['mageck_rra_count_min'],
        sorting_criteria_opt = config['mageck_rra_criteria'],
        test_threshold_opt = config['mageck_rra_pthreshold'],
        remove_from_opt = config['mageck_rra_remove'],
        control_sgrna_opt = lambda x: "--control-sgrna {controls_file}".format(controls_file=config['controls_file']) if config['controls_file'] != '' else ''
    log:
        "logs/{token}/MAGeCK/RRA/{treatment}_vs_{control}.log"
    shell:
        "mageck test \
            -k {input.count_table} \
            -t {params.treatment} \
            -c {params.control} \
            -n results/{wildcards.token}/MAGeCK_RRA/{wildcards.treatment}_vs_{wildcards.control}/{wildcards.treatment}_vs_{wildcards.control} \
            --gene-lfc-method {params.lfc_method_opt} \
            --adjust-method {params.ajd_method_opt} \
            --remove-zero-threshold {params.rm_zero_threshold_opt} \
            --sort-criteria {params.sorting_criteria_opt} \
            --gene-test-fdr-threshold {params.test_threshold_opt} \
            --remove-zero {params.remove_from_opt} \
            {params.control_sgrna_opt} &> {log}"
