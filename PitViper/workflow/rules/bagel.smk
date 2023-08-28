rule bagel_generate_count_matrix:
    """ Generate count matrix between two conditions for BAGEL2. """
    input:
        samples=config['tsv_file'],
        counts=rules.counts_filtering.output.normalized_filtered_counts#config['normalized_count_table']
    output:
        matrix="results/{token}/count_matrices/BAGEL2/{treatment}_vs_{control}_count_matrix.txt"

    log:
        "logs/{token}/BAGEL2/{treatment}_vs_{control}_count_matrix.log"
    message:
        "Generating count matrix for BAGEL2: {wildcards.treatment} vs {wildcards.control} \
        on files {input.samples} and {input.counts}. Output file: {output.matrix}"
    shell:
        "python3 workflow/scripts/bagel_count_files.py \
            --file {input.samples} \
            --counts {input.counts} \
            --directory results/{wildcards.token}/count_matrices/BAGEL2/ \
            --control {wildcards.control} \
            --treatment {wildcards.treatment} > {log}"


rule bagel_foldchange:
    "Run BAGEL2's script for foldchange calculation."
    input:
        count_table = rules.bagel_generate_count_matrix.output.matrix
    output:
        foldchange="results/{token}/BAGEL2/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL.foldchange"
    params:
        "results/{token}/BAGEL2/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL"
    log:
        "logs/{token}/BAGEL2/{treatment}_vs_{control}_foldchange.log"
    message:
        "Running BAGEL2's foldchange calculation for {wildcards.treatment} vs {wildcards.control} \
        on file {input.count_table}. Output file: {output.foldchange}"
    shell:
        "python workflow/scripts/bagel/BAGEL.py fc \
            -i {input.count_table} \
            -o {params} \
            -c 1 > {log}"


def bagel_bf_columns(wildcards):
    file = "results/{token}/BAGEL2/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL.foldchange".format(token=wildcards.token, treatment=wildcards.treatment, control=wildcards.control)
    content = pd.read_csv(file, sep="\t")
    out = ",".join([ str(i + 1) for i in range(len(content.columns)-2)])
    return out


rule bagel_bf:
    "Use BAGEL2's foldchanges for gene essentiality analysis."
    input:
        foldchange = rules.bagel_foldchange.output.foldchange
    output:
        bf = "results/{token}/BAGEL2/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.bf"
    params:
        nonessential = config['nonessentials'],
        essentials = config['essentials'],
        columns = bagel_bf_columns,
    log:
        "logs/{token}/BAGEL2/{treatment}_vs_{control}_bf.log"
    message:
        "Running BAGEL2's foldchange calculation for {wildcards.treatment} vs {wildcards.control} \
        on file {input.foldchange}. Output file: {output.bf}"
    shell:
        "python workflow/scripts/bagel/BAGEL.py bf \
            -i {input.foldchange} \
            -o {output.bf} \
            -e {params.essentials} \
            -n {params.nonessential} \
            -c {params.columns} > {log}"

rule bagel_pr:
    "Use BAGEL2's precision-recall."
    input:
        foldchange = rules.bagel_bf.output.bf
    output:
        pr = "results/{token}/BAGEL2/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.pr"
    params:
        nonessential = config['nonessentials'],
        essentials = config['essentials'],
    log:
        "logs/{token}/BAGEL2/{treatment}_vs_{control}_pr.log"
    message:
        "Running BAGEL2's foldchange calculation for {wildcards.treatment} vs {wildcards.control} \
        on file {input.foldchange}. Output file: {output.pr}"

    shell:
        "python workflow/scripts/bagel/BAGEL.py pr \
            -i {input.foldchange} \
            -o {output.pr} \
            -e {params.essentials} \
            -n {params.nonessential} > {log}"
