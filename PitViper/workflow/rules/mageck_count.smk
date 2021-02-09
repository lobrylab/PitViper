



rule mageck_count:
    input:
        fastqs = getFastqFiles,
        library = config['inputs']['library']
    output:
        file = "data/screen.count.txt"
    params:
        name = "data/screen",
        labels = getLabels,
        files = getFiles
    conda:
        "../envs/mageck.yaml"
    log:
        "logs/MAGeCK_counts.log"
    shell:
        "mageck count -l {input.library} -n {params.name} --sample-label {params.labels}  --fastq {params.files} &> {log}"


# rule mageck_mle:
#     """MAGeCK MLE implementation for gene essentiality analysis.
#         input: count table and design matrix
#         output: genes and sgrnas summaries
#         params: user must choose a normalization method"""
#     input:
#         designmat = rules.generate_design_matrix.output.matrix
#     output:
#         gene_summary = "results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt",
#         sgrna_summary = "results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/{treatment}_vs_{control}.sgrna_summary.txt"
#     params:
#         count_table = config['inputs']['count_table'],
#         name = "results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/{treatment}_vs_{control}",
#         method= config['MAGeCK']['MLE']['norm-method']
#     conda:
#         "../envs/mageck.yaml"
#     log:
#         "logs/{token}/MAGeCK/MLE/{treatment}_vs_{control}.log"
#     shell:
#         "mageck mle -k {params.count_table} -d {input.designmat} -n {params.name} --norm-method {params.method} &> {log}"
