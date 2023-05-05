

def getFastqMapping(wildcards):
    """Return name of fastq file."""
    samples_sheet = pd.read_csv(config['tsv_file'], sep="\t")
    fastq = samples_sheet.loc[samples_sheet["replicate"] == wildcards.sample].fastq
    return fastq

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
    log:
        "logs/" + config['token'] + "/mapping/bowtie_build.log"
    message:
        "Building bowtie index for {input}."
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
        fastq=getFastqMapping
    output:
        bam="results/{token}/bam/{sample}.bam"
    params:
        adap5=config['length_5_adapter'],
        adap3=config['length_3_adapter'],
        bowtie_mapping_method=config['bowtie_mapping_method'],
        index_base_name=index
    log:
        "logs/{token}/Bowtie_mapping/{sample}_bowtie_mapping.log"
    message:
        "Mapping {input.fastq} to {input.index}."
    shell:
        "bowtie2 \
            -x {params.index_base_name} \
            -U {input.fastq} \
            -5 {params.adap5} \
            -3 {params.adap3} \
            {params.bowtie_mapping_method} \
            --norc 2> {log} | samtools view -bS - > {output.bam}"
