

rule compute_mean:
    "Compute mean count of control condition for BAGEL."
    input:
        config['inputs']['count_table']
    output:
        count_table = 'data/{token}//{treatment}_vs_{control}_counts_for_BAGEL.txt'
    params:
        control=config['experimental_design']['control']
    shell:
        "python3 ./workflow/scripts/compute_mean.py {input} {params.control}"



rule bagel_foldchange:
    "Run BAGEL's script for foldchange calculation."
    input:
        count_table = rules.compute_mean.output.count_table
    output:
        foldchange="results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL.foldchange"
    params:
        "results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL"
    log:
        "logs/{token}/BAGEL/{treatment}_vs_{control}_foldchange.log"
    shell:
        "./workflow/scripts/bagel-for-knockout-screens-code/BAGEL-calc_foldchange.py \
            -i {input.count_table} \
            -o {params} \
            -c 10 > {log}"


rule bagel_bf:
    "Use BAGEL's foldchanges for gene essentiality analysis."
    input:
        foldchange = rules.bagel_foldchange.output.foldchange
    output:
        bf = "results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.bf"
    params:
        essentials = config['BAGEL']['essentials'],
        nonessential = config['BAGEL']['nonessential']
    log:
        "logs/{token}/BAGEL/{treatment}_vs_{control}_bf.log"
    shell:
        "./workflow/scripts/bagel-for-knockout-screens-code/BAGEL.py \
            -i {input.foldchange} \
            -o {output.bf} \
            -e {params.essentials} \
            -n {params.nonessential} \
            -c 7,8,9 > {log}"


rule bagel_essentials_genes:
    "Get essential genes estimated by BAGEL."
    input:
        rules.bagel_bf.output.bf
    output:
        "results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_bagel_essentials_genes.txt"
    shell:
        "gawk '(NR==1) {{print $0}} ($2 > 0 && NR!=1) {{print $0}}' {input} > {output}"


rule bagel_nonessentials_genes:
    "Get nonessential genes estimated by BAGEL."
    input:
        rules.bagel_bf.output.bf
    output:
        "results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_bagel_nonessentials_genes.txt"
    shell:
        "gawk '(NR==1) {{print $0}} ($2 < 0 && NR!=1) {{print $0}}' {input} > {output}"
