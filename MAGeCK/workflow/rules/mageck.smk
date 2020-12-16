

rule mageck_count_fast:
    "When starting from raw fastq files."
    input:
    output:
    params:
    conda:
        "../envs/mageck.yaml"
    log:
        "logs/{token}/mageck/count/{sample}"
    shell:
        "mageck count \
            -l library.txt \
            -n demo \
            --sample-label L1,CTRL  \
            --fastq test1.fastq test2.fastq"


rule mapping_bowtie2:
    input:

    output:

    params:

    conda:
        "../envs/bowtie2.yml"
    log:
        "logs/bowtie2/{sample}.log"
    threads: config['bowtie2']['threads']
    shell:
        "mageck count \
          -l {input.list} \
          -n {params.output_prefix} \
          --sample-label 'plasmid,ESC1' \
          --fastq ERR376998.bam ERR376999.bam"


rule mageck_RRA:
    "Simple treatment vs. control analysis."
    input:
    output:
    params:
    conda:
        "../envs/mageck.yaml"
    log:
        "logs/{token}/mageck/RRA/{sample}"
    shell:
        "mageck test \
            -k sample.txt \
            -t HL60.final,KBM7.final \
            -c HL60.initial,KBM7.initial \
              -n demo"


rule mageck_MLE:
    "Multiple sample comparison analysis."
    input:
    output:
    params:
    conda:
        "../envs/mageck.yaml"
    log:
        "logs/{token}/mageck/MLE/{sample}"
    shell:
        "mageck mle \
            -k leukemia.new.csv \
            -d designmat.txt \
            -n beta_leukemia"
