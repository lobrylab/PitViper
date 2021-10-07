

tsv = pd.read_csv(config['tsv_file'], sep="\t")
if len(tsv.columns) > 2:
    if "fastq" == tsv.columns[2]:
        ruleorder: mageck_count_fastq > mageck_count_bam
    elif "bam" == tsv.columns[2]:
        ruleorder: mageck_count_bam > mageck_count_fastq

rule mageck_count_fastq:
    input:
        fastqs = getFastqFiles,
        library = config['library_file']
    output:
        file = "results/" + config['token'] + "/" + config['count_table_file'],
        normalized = "results/" + config['token'] + "/" + config['normalized_count_table']
    params:
        name = "results/" + config['token'] + "/" + config['normalized_count_table'].split('.count_normalized.txt')[0],
        labels = getLabels,
        files = getFiles,
        revcom_opt = config['mageck_count_rev_comp'],
        normalization_opt = config['mageck_count_normalization'],
        length_opt = config['mageck_count_length'],
        align_all_opt = config['mageck_count_all_align'],
        count_N_opt = config['mageck_count_N']
    # conda:
    #     "../envs/mageck.yaml"
    log:
        "logs/mapping/MAGeCK_counts_fastq.log"
    shell:
        "mageck count -l {input.library} \
            -n {params.name} \
            --sample-label {params.labels} \
            --fastq {params.files} \
            --reverse-complement {params.revcom_opt} \
            --norm-method {params.normalization_opt} \
            --sgrna-len {params.length_opt} \
            --count-pair {params.align_all_opt} \
            --count-n {params.count_N_opt} &> {log}"


rule mageck_count_bam:
    input:
        bams = getBamFiles,
        library = config['library_file']
    output:
        file = "results/" + config['token'] + "/" + config['count_table_file'],
        normalized = "results/" + config['token'] + "/" + config['normalized_count_table']
    params:
        name = "results/" + config['token'] + "/" + config['normalized_count_table'].split('.count_normalized.txt')[0],
        labels = getLabels,
        files = getFiles,
        revcom_opt = config['mageck_count_rev_comp'],
        normalization_opt = config['mageck_count_normalization'],
        length_opt = config['mageck_count_length'],
        align_all_opt = config['mageck_count_all_align'],
        count_N_opt = config['mageck_count_N']
    # conda:
    #     "../envs/mageck.yaml"
    log:
        "logs/mapping/MAGeCK_counts_bam.log"
    shell:
        "mageck count -l {input.library} \
            -n {params.name} \
            --sample-label {params.labels} \
            --fastq {params.files} \
            --reverse-complement {params.revcom_opt} \
            --norm-method {params.normalization_opt} \
            --sgrna-len {params.length_opt} \
            --count-pair {params.align_all_opt} \
            --count-n {params.count_N_opt} &> {log}"
