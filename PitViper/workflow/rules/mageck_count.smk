

tsv = pd.read_csv(config['tsv_file'], sep="\t")
if len(tsv.columns) > 3:
    if config['bowtie_activate'] == "True":
        ruleorder: bowtie_mapping > mageck_count_fastq
        ruleorder: mageck_count_bam > mageck_count_fastq
        ruleorder: mageck_count_bam > MAGeCK_counts_normalize
        ruleorder: MAGeCK_counts_normalize > mageck_count_fastq
    elif config['mageck_count_activate'] == "True":
        if "fastq" in tsv.columns:
            ruleorder: mageck_count_fastq > mageck_count_bam
            ruleorder: mageck_count_fastq > MAGeCK_counts_normalize
        elif "bam" in tsv.columns:
            ruleorder: mageck_count_bam > mageck_count_fastq
            ruleorder: mageck_count_bam > MAGeCK_counts_normalize
else:
    ruleorder: MAGeCK_counts_normalize > mageck_count_bam
    ruleorder: MAGeCK_counts_normalize > mageck_count_fastq


rule mageck_count_fastq:
    input:
        fastqs = getFastqFiles,
        library = config['library_file']
    output:
        file = config['count_table_file'],
        normalized = config['normalized_count_table']
    params:
        name = config['count_table_file'].split('.count.txt')[0],
        labels = getLabels,
        files = getFiles,
        revcom_opt = lambda x: "--reverse-complement" if config['mageck_count_rev_comp'] == 'True' else '',
        normalization_opt = config['mageck_count_normalization'],
        length_opt = config['mageck_count_length'],
        align_all_opt = lambda x: "--count-pair True" if config['mageck_count_all_align'] == 'True' else '',
        count_N_opt = lambda x: "--count-n" if config['mageck_count_N'] == 'True' else ''
    log:
        "logs/mapping/MAGeCK_counts_fastq.log"
    shell:
        "mageck count -l {input.library} \
            -n {params.name} \
            --sample-label {params.labels} \
            --fastq {params.files} \
            {params.revcom_opt} \
            --norm-method {params.normalization_opt} \
            --sgrna-len {params.length_opt} \
            {params.align_all_opt} \
            {params.count_N_opt} &> {log}"


rule mageck_count_bam:
    input:
        bams = getBamFiles,
        library = config['library_file']
    output:
        file = config['count_table_file'],
        normalized = config['normalized_count_table']
    params:
        name = config['count_table_file'].split('.count.txt')[0],
        labels = getLabels,
        files = getFiles,
        revcom_opt = lambda x: "--reverse-complement" if config['mageck_count_rev_comp'] == 'True' else '',
        normalization_opt = config['mageck_count_normalization'],
        length_opt = config['mageck_count_length'],
        align_all_opt = lambda x: "True" if config['mageck_count_all_align'] == 'True' else 'False',
        count_N_opt = lambda x: "--count-n" if config['mageck_count_N'] == 'True' else ''
    log:
        "logs/mapping/MAGeCK_counts_bam.log"
    shell:
        "mageck count -l {input.library} \
            -n {params.name} \
            --sample-label {params.labels} \
            --fastq {params.files} \
            {params.revcom_opt} \
            --norm-method {params.normalization_opt} \
            --sgrna-len {params.length_opt} \
            --count-pair {params.align_all_opt} \
            {params.count_N_opt} &> {log}"
