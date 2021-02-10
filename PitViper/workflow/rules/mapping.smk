



rule library_convert:
    input:
        csv=config['inputs']['library']
    output:
        fa=config['inputs']['library'] + ".fa"
    shell:
        "./workflow/scripts/convertLibrary.sh {input.csv} {output.fa}"


rule bowtie_build:
    input:
        rules.library_convert.output.fa
    output:
        index=config['inputs']['library'][0:-4] + ".1.bt2"
    params:
        config['inputs']['library'][0:-4]
    conda:
        "../envs/bowtie.yaml"
    log:
        "logs/mapping/bowtie_build.log"
    shell:
        "bowtie2-build {input} {params} > {log}"


rule bowtie_mapping:
    input:
        index=rules.bowtie_build.output.index,
        fastq="path/to/fastq/{sample}.fastq"
    output:
        bam="path/to/bam/{sample}.bam"
    params:
        adap5="",
        adap3=""
    conda:
        "../envs/bowtie.yaml"
    log: ""
    shell:
        "bowtie2 \
            -x {input.index} \
            -U {input.fastq} \
            -5 {params.adap5} \
            -3 {params.adap3} \
            --norc | samtools view -bS - > {output.bam}"
