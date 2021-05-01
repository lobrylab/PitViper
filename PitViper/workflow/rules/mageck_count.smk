

tsv = pd.read_csv(config['inputs']['tsv'], sep="\t")
if "fastq" == tsv.columns[2]:
    ruleorder: mageck_count_fastq > mageck_count_bam
elif "bam" == tsv.columns[2]:
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
