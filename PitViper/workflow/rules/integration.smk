

def generatedResults(wildcards):
    results = []
    comparaisons = get_all_pairwise_comparisons()
    token = config['token']
    for comparaison in comparaisons:
        if (config['bagel_activate'] == 'True'):
            results.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.bf".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
            results.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_bagel_essentials_genes.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
            results.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_bagel_essentials_genes.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
            results.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_bagel_nonessentials_genes.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['mageck_mle_activate'] == 'True'):
            results.append("results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['mageck_rra_activate'] == 'True'):
            results.append("results/{token}/MAGeCK_RRA/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['crisphiermix_activate'] == 'True'):
            results.append("results/{token}/CRISPhieRmix/{treatment}_vs_{control}/{treatment}_vs_{control}.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['gsea_activate'] == 'True'):
            results.append("results/{token}/GSEA-like/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_GSEA-like.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['filtering_activate'] == 'True'):
            results.append("results/{token}/in_house_method/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_in-house.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
    return results


rule genes_integration:
    input:
        generatedResults
    output:
        "results/" + config['token'] + "/ExecutionComplete.txt"
    params:
        config['token']
    log:
        notebook="results/" + config['token'] + "/Report.ipynb"
    notebook:
        "../notebooks/Report_template.py.ipynb"
