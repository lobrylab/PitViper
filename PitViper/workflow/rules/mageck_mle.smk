rule generate_design_matrix:
    """Create design matrix for treatment versus control."""
    input:
        config['tsv_file']
    output:
        matrix="results/{token}/design_matrices/MAGeCK/{treatment}_vs_{control}_design_matrix.txt"
    log:
        "logs/{token}/MAGeCK/MLE/{treatment}_vs_{control}_design_matrix.log"
    message:
        "Generating design matrix for MAGeCK MLE analysis. From {input} to {output.matrix}."
    shell:
        "python3 workflow/scripts/readSamples.py --file {input} --directory results/{wildcards.token}/design_matrices/MAGeCK/ --control {wildcards.control} --treatment {wildcards.treatment}"


rule mageck_mle:
    """MAGeCK MLE implementation for gene essentiality analysis.
        input: count table and design matrix
        output: genes and sgrnas summaries
        params: user must choose a normalization method"""
    input:
        designmat = rules.generate_design_matrix.output.matrix,
        count_table = rules.counts_filtering.output.raw_filtered_counts
    output:
        gene_summary = "results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt",
        sgrna_summary = "results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/{treatment}_vs_{control}.sgrna_summary.txt"
    params:
        name = "results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/{treatment}_vs_{control}",
        method= config['mageck_mle_normalization'],
        adj_opt = config['mageck_mle_adj'],
        mean_var_opt = config['mageck_mle_mean_var'],
        outliers_opt = lambda x: "--remove-outliers" if config['mageck_mle_outliers'] == 'True' else '',
        perm_round_opt = config['mageck_mle_perm_N'],
        no_perm_group_opt = lambda x: "--no-permutation-by-group" if config['mageck_mle_perm_all'] == 'True' else '',
        control_sgrna_opt = lambda x: "--control-sgrna {controls_file}".format(controls_file=config['controls_file']) if config['controls_file'] != '' else ''
    log:
        "logs/{token}/MAGeCK/MLE/{treatment}_vs_{control}.log"
    message:
        "Running MAGeCK MLE for {params.name}. From {input.count_table} and {input.designmat} to {output.gene_summary} and {output.sgrna_summary}."
    shell:
        "mageck mle -k {input.count_table} \
            -d {input.designmat} \
            -n {params.name} \
            --threads 10 \
            --adjust-method {params.adj_opt} \
            --genes-varmodeling {params.mean_var_opt} \
            --norm-method {params.method} \
            {params.outliers_opt} \
            --permutation-round {params.perm_round_opt} \
            {params.no_perm_group_opt} \
            {params.control_sgrna_opt} &> {log}"
