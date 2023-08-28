

def generatedResults(wildcards):
    results = []
    comparaisons = get_all_pairwise_comparisons()
    token = config['token']
    for comparaison in comparaisons:
        if (config['bagel_activate'] == 'True'):
            results.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.bf".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
            results.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.pr".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['mageck_mle_activate'] == 'True'):
            results.append("results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['mageck_rra_activate'] == 'True'):
            results.append("results/{token}/MAGeCK_RRA/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['crisphiermix_activate'] == 'True'):
            results.append("results/{token}/CRISPhieRmix/{treatment}_vs_{control}/{treatment}_vs_{control}.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['ssrea_activate'] == 'True'):
            results.append("results/{token}/SSREA/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_SSREA.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['directional_scoring_method_activate'] == 'True'):
            results.append("results/{token}/directional_scoring_method/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_directional_scoring_method.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        results.append(config['normalized_count_table'])
    return results


rule genes_integration:
    input:
        generatedResults
    output:
        "results/" + config['token'] + "/Report.ipynb"
    params:
        config['token']
    log:
        notebook="results/" + config['token'] + "/Report.ipynb"
    notebook:
        "../notebooks/Report_template.py.ipynb"
