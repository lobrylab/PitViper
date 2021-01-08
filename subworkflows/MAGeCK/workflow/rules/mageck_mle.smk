# rule rule:
#     input:
#     output:
#     params:
#     conda:
#     log:
#     notebook:
#     shell:


rule mageck_mle:
    """MAGeCK MLE """
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
        "../envs/mageck.yml"
    log:
        logs="logs/{token}/MAGeCK/MLE/{name}.log"
    shell:
        "mageck mle -k {input.count_table} -d {input.designmat} -n {params.name} --norm-method {params.method} 2>1 {log}"


rule notebooks:
    input:
        gene_summary=rules.mageck_mle.output.gene_summary,
        sgrna_summary=rules.mageck_mle.output.sgrna_summary
    output:
        "results/{token}/notebooks/MAGeCK_MLE_{name}.txt"
    conda:
        "../envs/jupyter.yml"
    log:
        # optional path to the processed notebook
        notebook="logs/notebooks/{token}/{name}_processed_notebook.ipynb"
    notebook:
        "notebooks/{token}/{name}.py.ipynb"
