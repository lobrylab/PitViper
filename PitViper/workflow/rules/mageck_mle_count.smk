



rule mageck_mle:
    """MAGeCK MLE implementation for gene essentiality analysis.
        input: count table and design matrix
        output: genes and sgrnas summaries
        params: user must choose a normalization method"""
    input:
        count_table = config['inputs']['count_table'],
        designmat = config['MAGeCK']['designmat']
    output:
        gene_summary = "results/{token}/MAGeCK_MLE/{name}.gene_summary.txt",
        sgrna_summary = "results/{token}/MAGeCK_MLE/{name}.sgrna_summary.txt"
    params:
        name = "results/{token}/MAGeCK_MLE/{name}",
        method= config['MAGeCK']['MLE']['norm-method']
    conda:
        "../envs/mageck.yaml"
    log:
        logs="logs/{token}/MAGeCK/MLE/{name}.log"
    shell:
        "mageck mle -k {input.count_table} -d {input.designmat} -n {params.name} --norm-method {params.method} > {log}"


rule mageck_mle_notebooks:
    input:
        gene_summary=rules.mageck_mle.output.gene_summary,
        sgrna_summary=rules.mageck_mle.output.sgrna_summary
    output:
        txt="results/{token}/MAGeCK_MLE/MAGeCK_MLE_{name}.txt",
        pdf="results/{token}/MAGeCK_MLE/MAGeCK_MLE_{name}_{condition}_genes_beta_plot.html"
    conda:
        "../envs/jupyter.yaml"
    log:
        notebook="notebooks/{token}/{name}_mageck_mle_processed_notebook.ipynb"
    notebook:
        "../notebooks/MAGeCK_MLE.py.ipynb"
