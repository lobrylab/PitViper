



rule mapping_bowtie2:
    input:
        mate1="",
        mate2=""
    output:
        bam="",
        stat=""
    params:
        index="",
        mismatches_allowed=config['bowtie2']['mismatches_allowed']
    conda:
        "../envs/bowtie2.yml"
    log:
        "logs/bowtie2/{sample}.log"
    threads: config['bowtie2']['threads']
    shell:
        "bowtie2 \
          -x {params.index} \
          -1 {input.mate1} \
          -2 {input.mate2} \
          --very-sensitive \
          -N {params.mismatches_allowed} \
          --end-to-end \
          -p {threads} 2> {output.stat} \
          | samtools view -S -b - | samtools sort -@ {threads} - > {output.bam}"

"bowtie2 \
  -x {params.index} \
  -U {input.fastq} \
  -5 23 \
  -3 8 \
  --norc 2> {output.stat}\
  | samtools view -bS - > {output.sam}"
