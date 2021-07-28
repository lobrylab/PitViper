



rule crisphiermix_generate_count_matrix:
    """ Generate a count matrix between two conditions for CRISPhieRmix. """
    input:
        samples=config['inputs']['tsv'],
        counts="results/{token}/" + config['inputs']['count_table']
    output:
        matrix="results/{token}/count_matrices/CRISPhieRmix/{treatment}_vs_{control}_count_matrix.txt"
    conda:
        "../envs/commons.yaml"
    log:
        "logs/{token}/CRISPhieRmix/{treatment}_vs_{control}_count_matrix.log"
    shell:
        "python3 workflow/scripts/crisphiermix_count_files.py \
            --file {input.samples} \
            --counts {input.counts} \
            --directory results/{wildcards.token}/count_matrices/CRISPhieRmix/ \
            --control {wildcards.control} \
            --treatment {wildcards.treatment} > {log}"


rule CRISPhieRmix:
    """
        CRISPhieRmix implementation for gene essentiality analysis.
        input: DESeq2 results
        output: genes summaries
    """
    input:
        deseq2_table=rules.DESeq2_counts.output.deseq2_out,
        count_table="results/{token}/" + config['inputs']['count_table']
    output:
        gene_summary="results/{token}/CRISPhieRmix/{treatment}_vs_{control}/{treatment}_vs_{control}.txt",
    params:
        n_control=getControlIdsLen,
        n_treatment=getTreatmentIdsLen
    conda:
        "../envs/crisphiermix.yaml"
    log:
        "logs/{token}/CRISPhieRmix/{treatment}_vs_{control}.log"
    benchmark:
        "benchmarks/{token}/CRISPhieRmix/{treatment}_vs_{control}.tsv"
    script:
        "../../workflow/scripts/run_crisphiermix.R"
