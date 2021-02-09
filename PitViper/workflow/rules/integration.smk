



rule genes_integration:
    input:
        mageck_mle=rules.mageck_mle.output.gene_summary,
        mageck_rra=rules.mageck_rra.output.gene_summary,
        bagel=rules.bagel_bf.output.bf,
        crisphiermix=rules.CRISPhieRmix.output.gene_summary
    output:
        mageck_mle_genes="results/{token}/Integration/{treatment}_vs_{control}/Genes/MAGeck_MLE_{treatment}_vs_{control}.txt",
        mageck_rra_genes="results/{token}/Integration/{treatment}_vs_{control}/Genes/MAGeck_RRA_{treatment}_vs_{control}.txt",
        bagel_genes="results/{token}/Integration/{treatment}_vs_{control}/Genes/BAGEL_{treatment}_vs_{control}.txt",
        crisphiermix="results/{token}/Integration/{treatment}_vs_{control}/Genes/CRISPhieRmix_{treatment}_vs_{control}.txt"
    conda:
        "../envs/jupyter.yaml"
    log:
        notebook="notebooks/{token}/{treatment}_vs_{control}_genes_integration_notebook.ipynb"
    notebook:
        "../notebooks/Genes_Integration.py.ipynb"
