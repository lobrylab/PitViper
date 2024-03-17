rule bagel_generate_count_matrix:
    """ Generate count matrix between two conditions for BAGEL. """
    input:
        samples=config['tsv_file'],
        counts=rules.counts_filtering.output.normalized_filtered_counts
    output:
        matrix="results/{token}/count_matrices/BAGEL/{treatment}_vs_{control}_count_matrix.txt"

    log:
        "logs/{token}/BAGEL/{treatment}_vs_{control}_count_matrix.log"
    message:
        "Generating count matrix for BAGEL: {wildcards.treatment} vs {wildcards.control} \
        on files {input.samples} and {input.counts}. Output file: {output.matrix}"
    shell:
        "python3 workflow/scripts/bagel_count_files.py \
            --file {input.samples} \
            --counts {input.counts} \
            --directory results/{wildcards.token}/count_matrices/BAGEL/ \
            --control {wildcards.control} \
            --treatment {wildcards.treatment} > {log}"


rule bagel_foldchange:
    "Run BAGEL's script for foldchange calculation."
    input:
        count_table = rules.bagel_generate_count_matrix.output.matrix
        # count_table = rules.counts_filtering.output.normalized_filtered_counts
    output:
        foldchange="results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL.foldchange"
    params:
        out_dir="results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL",
        control_cols=1,  #getControlIds,
    log:
        "logs/{token}/BAGEL/{treatment}_vs_{control}_foldchange.log"
    message:
        "Running BAGEL's foldchange calculation for {wildcards.treatment} vs {wildcards.control} \
        on file {input.count_table}. Output file: {output.foldchange}"
    shell:
        "python workflow/scripts/bagel/BAGEL.py fc \
            -i {input.count_table} \
            -o {params.out_dir} \
            -c {params.control_cols} > {log}"


rule bagel_bf:
    "Use BAGEL's foldchanges for gene essentiality analysis."
    input:
        foldchange = rules.bagel_foldchange.output.foldchange
    output:
        bf = "results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.bf"
    params:
        nonessential = config['nonessentials'],
        essentials = config['essentials'],
        columns = getTreatmentIds,
    log:
        "logs/{token}/BAGEL/{treatment}_vs_{control}_bf.log"
    message:
        "Running BAGEL's foldchange calculation for {wildcards.treatment} vs {wildcards.control} \
        on file {input.foldchange}. Output file: {output.bf}"
    shell:
        "python workflow/scripts/bagel/BAGEL.py bf \
            -i {input.foldchange} \
            -o {output.bf} \
            -e {params.essentials} \
            -n {params.nonessential} \
            -s 123 \
            -c {params.columns} > {log}"

rule bagel_pr:
    "Use BAGEL's precision-recall."
    input:
        foldchange = rules.bagel_bf.output.bf
    output:
        pr = "results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.pr"
    params:
        nonessential = config['nonessentials'],
        essentials = config['essentials'],
    log:
        "logs/{token}/BAGEL/{treatment}_vs_{control}_pr.log"
    message:
        "Running BAGEL's foldchange calculation for {wildcards.treatment} vs {wildcards.control} \
        on file {input.foldchange}. Output file: {output.pr}"

    shell:
        "python workflow/scripts/bagel/BAGEL.py pr \
            -i {input.foldchange} \
            -o {output.pr} \
            -e {params.essentials} \
            -n {params.nonessential} > {log}"
