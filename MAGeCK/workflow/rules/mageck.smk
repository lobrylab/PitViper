

rule mageck_count:
    input:
    output:
    params:
    conda:
        "../envs/mageck.yaml"
    logs:
        "logs/{token}/mageck/count/{sample}"
    shell:
        "mageck count \
            -l library.txt \
            -n {params} \
            --sample-label L1,CTRL  \
            --fastq test1.fastq test2.fastq"


rule mageck_RRA:
    input:
        "results/{token}/counts/"
    output:
    params:
    conda:
        "../envs/mageck.yaml"
    logs:
        "logs/{token}/mageck/RRA/{sample}"
    shell:
        "mageck test \
            -k sample.txt \
            -t HL60.final,KBM7.final \
            -c HL60.initial,KBM7.initial \
              -n demo"


rule mageck_MLE:
    input:
    output:
    params:
    conda:
        "../envs/mageck.yaml"
    logs:
        "logs/{token}/mageck/MLE/{sample}"
    shell:
        "mageck mle \
            -k leukemia.new.csv \
            -d designmat.txt \
            -n beta_leukemia"
