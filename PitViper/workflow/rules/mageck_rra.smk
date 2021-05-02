


rule mageck_rra:
    """MAGeCK RRA implementation for gene essentiality analysis.
        input: count table and design matrix
        output: genes and sgrnas summaries
        params: user must choose a normalization method"""
    input:
        count_table="results/{token}/" + config['inputs']['normalized_count_table']
    output:
        gene_summary = "results/{token}/MAGeCK_RRA/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt",
    params:
        name = "results/{token}/MAGeCK_RRA/{treatment}_vs_{control}/{treatment}_vs_{control}",
        method = config['MAGeCK']['MLE']['norm-method'],
        treatment = getTreatmentIds,
        control = getControlIds
    conda:
        "../envs/mageck.yaml"
    log:
        "logs/{token}/MAGeCK/RRA/{treatment}_vs_{control}.log"
    shell:
        "mageck test \
            -k {input.count_table} \
            -t {params.treatment} \
            -c {params.control} \
            -n results/{wildcards.token}/MAGeCK_RRA/{wildcards.treatment}_vs_{wildcards.control}/{wildcards.treatment}_vs_{wildcards.control} &> {log}"
