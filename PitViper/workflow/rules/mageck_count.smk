

if config['start_from'] == 'fastq':
    ruleorder: mageck_count_fastq > mageck_count_bam
elif config['start_from'] == 'bam':
    ruleorder: mageck_count_bam > mageck_count_fastq

rule mageck_count_fastq:
    input:
        fastqs = getFastqFiles,
        library = config['inputs']['library']
    output:
        file = config['inputs']['count_table']
    params:
        name = config['inputs']['count_table'].split('.count.txt')[0],
        labels = getLabels,
        files = getFiles
    conda:
        "../envs/mageck.yaml"
    log:
        "logs/mapping/MAGeCK_counts_fastq.log"
    shell:
        "mageck count -l {input.library} -n {params.name} --sample-label {params.labels}  --fastq {params.files} &> {log}"


rule mageck_count_bam:
    input:
        bams = getBamFiles,
        library = config['inputs']['library']
    output:
        file = config['inputs']['count_table']
    params:
        name = config['inputs']['count_table'].split('/')[-1].split('.count.txt')[0],
        labels = getLabels,
        files = getFiles
    conda:
        "../envs/mageck.yaml"
    log:
        "logs/mapping/MAGeCK_counts_bam.log"
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
