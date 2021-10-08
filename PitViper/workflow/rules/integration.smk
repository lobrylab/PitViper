



# rule genes_integration:
#     input:
#         mageck_mle=rules.mageck_mle.output.gene_summary,
#         mageck_rra=rules.mageck_rra.output.gene_summary,
#         bagel=rules.bagel_bf.output.bf,
#         crisphiermix=rules.CRISPhieRmix.output.gene_summary
#     output:
#         mageck_mle_genes="results/{token}/Integration/{treatment}_vs_{control}/Genes/MAGeck_MLE_{treatment}_vs_{control}.txt",
#         mageck_rra_genes="results/{token}/Integration/{treatment}_vs_{control}/Genes/MAGeck_RRA_{treatment}_vs_{control}.txt",
#         bagel_genes="results/{token}/Integration/{treatment}_vs_{control}/Genes/BAGEL_{treatment}_vs_{control}.txt",
#         crisphiermix="results/{token}/Integration/{treatment}_vs_{control}/Genes/CRISPhieRmix_{treatment}_vs_{control}.txt"
#     conda:
#         "../envs/jupyter.yaml"
#     log:
#         notebook="notebooks/{token}/{treatment}_vs_{control}_genes_integration_notebook.ipynb"
#     notebook:
#         "../notebooks/Genes_Integration.py.ipynb"


def generatedResults(wildcards):
    results = []
    comparaisons = get_all_pairwise_comparaisons()

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


# rule integration:
#     input:
#         generatedResults
#     output:
#         # nb="results/" + config['token'] + "/reports/PitViper_report_" + config['token'] + ".ipynb"
#         "test.txt"
#     log:
#         notebook="logs/notebooks/processed_notebook.ipynb"
#         # notebook="results/" + config['token'] + "/reports/PitViper_report_" + config['token'] + ".ipynb"
#     notebook:
#         "notebooks/hello.py.ipynb"
#         # "workflow/notebooks/Rapport_PitViper_copie.ipynb"


rule genes_integration:
    input:
        generatedResults
    output:
        "results/" + config['token'] + "/ExecutionComplete.txt"
    params:
        config['token']
    # conda:
    #     "../envs/jupyter.yaml"
    log:
        notebook="results/" + config['token'] + "/Report.ipynb"
    notebook:
        "../notebooks/Report_template.py.ipynb"


# rule generate_report:
#     input:
#         generatedResults
#     output:
#         notebook="results/" + config['token'] + "/reports/PitViper_report_" + config['token'] + ".ipynb"
#     params:
#         template="workflow/notebooks/Rapport_PitViper.ipynb",
#         token=config['token']
#     # conda:
#     #     "../envs/jupyter.yaml"
#     log:
#         "logs/" + config['token'] + "/reports/PitViper_report.log"
#     shell:
#         "papermill {params.template} {output.notebook} \
#             -p token {params.token}"


# rule cp_report:
#     input:
#         generatedResults,
#         template="workflow/notebooks/Rapport_PitViper.ipynb"
#     output:
#         notebook="results/" + config['token'] + "/reports/PitViper_report_" + config['token'] + ".ipynb"
#     log:
#         "logs/" + config['token'] + "/reports/PitViper_report.log"
#     shell:
#         "cp {input.template} {output.notebook}"

# rule visualization:
#     input:
#         generatedResults
#     output:
#         nb="results/" + config['token'] + "/reports/integration.ipynb"
#     params:
#         in_nb="workflow/notebooks/Reports_Integration.py.ipynb",
#         token=config['token']
#     conda:
#         "../envs/jupyter.yaml"
#     log:
#         "logs/" + config['token'] + "/Integration/integration_test.log"
#     shell:
#         "papermill {params.in_nb} {output.nb} -p mageck_mle_outputs results/{params.token}/MAGeCK_MLE/ -p mageck_rra_outputs results/{params.token}/MAGeCK_RRA/"
