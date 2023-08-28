



rule crisphiermix_generate_count_matrix:
    """ Generate a count matrix between two conditions for CRISPhieRmix. """
    input:
        samples=config['tsv_file'],
        counts=rules.counts_filtering.output.normalized_filtered_counts#config['normalized_count_table']
    output:
        matrix="results/{token}/count_matrices/CRISPhieRmix/{treatment}_vs_{control}_count_matrix.txt"
    log:
        "logs/{token}/CRISPhieRmix/{treatment}_vs_{control}_count_matrix.log"
    message:
        "Generating count matrix for CRISPhieRmix: {wildcards.treatment} vs {wildcards.control}, using {input.counts} and {input.samples}."
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
        count_table=rules.counts_filtering.output.normalized_filtered_counts#config['normalized_count_table']
    output:
        gene_summary="results/{token}/CRISPhieRmix/{treatment}_vs_{control}/{treatment}_vs_{control}.txt",
    params:
        n_control=getControlIdsLen,
        n_treatment=getTreatmentIdsLen,
        control_sgrna_opt = lambda x: config['controls_file'] if config['controls_file'] != '' else '',
        screen_type = config['crisphiermix_type'],
        mu = config['crisphiermix_mu'],
        bimodal = config['crisphiermix_bimodal'],
    log:
        "logs/{token}/CRISPhieRmix/{treatment}_vs_{control}.log"
    benchmark:
        "benchmarks/{token}/CRISPhieRmix/{treatment}_vs_{control}.tsv"
    message:
        "Running CRISPhieRmix for {wildcards.treatment} vs {wildcards.control}."
    script:
        "../../workflow/scripts/run_crisphiermix.R"
