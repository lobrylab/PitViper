



rule generate_design_matrix:
    """Create design matrix for treatment versus control."""
    input:
        config['tsv_file']
    output:
        matrix="results/{token}/design_matrices/MAGeCK/{treatment}_vs_{control}_design_matrix.txt"
    # conda:
    #     "../envs/commons.yaml"
    log:
        "logs/{token}/MAGeCK/MLE/{treatment}_vs_{control}_design_matrix.log"
    shell:
        "python3 workflow/scripts/readSamples.py --file {input} --directory results/{wildcards.token}/design_matrices/MAGeCK/ --control {wildcards.control} --treatment {wildcards.treatment}"


rule mageck_mle:
    """MAGeCK MLE implementation for gene essentiality analysis.
        input: count table and design matrix
        output: genes and sgrnas summaries
        params: user must choose a normalization method"""
    input:
        designmat = rules.generate_design_matrix.output.matrix,
        count_table = config['count_table_file']
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
        no_perm_group_opt = lambda x: "--no-permutation-by-group" if config['mageck_mle_perm_all'] == 'True' else ''
    # conda:
    #     "../envs/mageck.yaml"
    log:
        "logs/{token}/MAGeCK/MLE/{treatment}_vs_{control}.log"
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
            {params.no_perm_group_opt} &> {log}"


# rule mageck_mle_notebooks:
#     """ Generate a jupyter notebook for data analysis of MAGeCK MLE results. """
#     input:
#         gene_summary=rules.mageck_mle.output.gene_summary,
#         sgrna_summary=rules.mageck_mle.output.sgrna_summary
#     output:
#         txt="results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/MAGeCK_MLE_{treatment}_vs_{control}.txt"
#     conda:
#         "../envs/jupyter.yaml"
#     log:
#         notebook="notebooks/{token}/{treatment}_vs_{control}_mageck_mle_processed_notebook.ipynb"
#     notebook:
#         "../notebooks/MAGeCK_MLE.py.ipynb"
