



rule library_convert:
    input:
        csv=config['library_file']
    output:
        fa=config['library_file'] + ".fa"
    shell:
        "./workflow/scripts/convertLibrary.sh {input.csv} {output.fa}"


rule bowtie_build:
    input:
        rules.library_convert.output.fa
    output:
        index=config['library_file'][0:-4] + ".1.bt2"
    params:
        config['library_file'][0:-4]
    # conda:
    #     "../envs/bowtie.yaml"
    log:
        "logs/mapping/bowtie_build.log"
    shell:
        "bowtie2-build {input} {params} > {log}"


def index(wildcards):
    full_name = config['library_file']
    m = re.match("(^.+)\.\w+$", full_name)
    if m:
        return m.group(1)


rule bowtie_mapping:
    input:
        index=rules.bowtie_build.output.index,
        fastq="../data/fastq/{sample}.fastq.gz"
    output:
        bam="../data/bam/{sample}.bam"
    params:
        adap5="20",
        adap3="10",
        index_base_name=index
    # conda:
    #     "../envs/bowtie.yaml"
    log:
        "logs/Bowtie_mapping/{sample}_bowtie_mapping.log"
    shell:
        "bowtie2 \
            -x {params.index_base_name} \
            -U {input.fastq} \
            -5 {params.adap5} \
            -3 {params.adap3} \
            --norc | samtools view -bS - > {output.bam}"
