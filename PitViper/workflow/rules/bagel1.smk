rule BAGEL1_foldchange:
    "Run BAGEL's script for foldchange calculation."
    input:
        count_table = "results/{token}/count_matrices/BAGEL/{treatment}_vs_{control}_count_matrix.txt"
    output:
        foldchange="results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL1.foldchange"
    params:
        "results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL1"
    conda:
        "../envs/bagel.yaml"
    log:
        "logs/{token}/BAGEL1/{treatment}_vs_{control}_foldchange.log"
    message:
        "Running BAGEL's foldchange calculation for {wildcards.treatment} vs {wildcards.control} \
        on file {input.count_table}. Output file: {output.foldchange}"
    shell:
        "python workflow/scripts/bagel-for-knockout-screens-code/BAGEL-calc_foldchange.py \
            -i {input.count_table} \
            -o {params} \
            -c 1 > {log}"


def bagel_bf_columns(wildcards):
    file = "results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL1.foldchange".format(token=wildcards.token, treatment=wildcards.treatment, control=wildcards.control)
    content = pd.read_csv(file, sep="\t")
    out = ",".join([ str(i + 1) for i in range(len(content.columns)-2)])
    return out


rule BAGEL1_bf:
    "Use BAGEL's foldchanges for gene essentiality analysis."
    input:
        foldchange = rules.BAGEL1_foldchange.output.foldchange
    output:
        bf = "results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL1_output.bf"
    params:
        nonessential = config['nonessentials'],
        essentials = config['essentials'],
        columns = bagel_bf_columns,
    log:
        "logs/{token}/BAGEL1/{treatment}_vs_{control}_bf.log"
    conda:
        "../envs/bagel.yaml"
    message:
        "Running BAGEL's foldchange calculation for {wildcards.treatment} vs {wildcards.control} \
        on file {input.foldchange}. Output file: {output.bf}"
    shell:
        "python workflow/scripts/bagel-for-knockout-screens-code/BAGEL.py \
            -i {input.foldchange} \
            -o {output.bf} \
            -e {params.essentials} \
            -n {params.nonessential} \
            -c {params.columns} > {log}"