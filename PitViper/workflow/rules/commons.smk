
def getFastqFiles(wildcards):
    """Return a list of experiment's files  in fastq format."""
    samples_sheet = pd.read_csv(config['tsv_file'], sep="\t")
    fastqs = samples_sheet.fastq.values
    return fastqs


def getBamFiles(wildcards):
    """Return a list of experiment's files in bam format."""
    samples_sheet = pd.read_csv(config['tsv_file'], sep="\t")
    bams = samples_sheet.bam.values
    return bams


def getLabels(wildcards):
    """Return concatenation of experiment's labels. Needed for MAGeCK count software."""
    samples_sheet = pd.read_csv(config['tsv_file'], sep="\t")
    labels = []
    for i in range(0, len(samples_sheet.index)):
        labels.append(samples_sheet.loc[i].replicate)
    return ",".join(labels)


def getFiles(wildcards):
    """Return concatenation of experiment's files. Needed for MAGeCK count software."""
    samples_sheet = pd.read_csv(config['tsv_file'], sep="\t")
    if 'fastq' in samples_sheet.columns:
        fastqs = []
        for i in range(0, len(samples_sheet.index)):
            fastqs.append(samples_sheet.loc[i].fastq)
        return " ".join(fastqs)
    elif 'bam' in samples_sheet.columns:
        bams = []
        for i in range(0, len(samples_sheet.index)):
            bams.append(samples_sheet.loc[i].bams)
        return " ".join(bams)


def getTreatmentIdsLen(wildcards):
    samples_file = pd.read_csv(config['tsv_file'], sep="\t")
    list_ids = samples_file.loc[samples_file["condition"] == wildcards.treatment].values
    return len(list_ids)


def getControlIdsLen(wildcards):
    samples_file = pd.read_csv(config['tsv_file'], sep="\t")
    list_ids = samples_file.loc[samples_file["condition"] == wildcards.control].values
    return len(list_ids)


def getTreatmentIds(wildcards):
    samples_file = pd.read_csv(config['tsv_file'], sep="\t")
    list_ids = samples_file.loc[samples_file.condition == wildcards.treatment].replicate.values
    return ",".join(list_ids)


def getControlIds(wildcards):
    samples_file = pd.read_csv(config['tsv_file'], sep="\t")
    list_ids = samples_file.loc[samples_file.condition == wildcards.control].replicate.values
    return ",".join(list_ids)


def get_all_pairwise_comparisons():
    design = pd.read_csv(config['tsv_file'], sep="\t")
    design_dict = {}
    comparisons = []
    for (k, v) in zip(design.order, design.condition):
        if not k in design_dict.keys():
            design_dict[k] = []
        if not v in design_dict[k]:
            design_dict[k].append(v)
    for i in range(len(design_dict)-1):
        for control in design_dict[i]:
            for k in range(i+1, len(design_dict)):
                for treatment in design_dict[k]:
                    comparisons.append({'treatment': treatment, 'control': control})
    return comparisons


def get_pipeline_outputs(wildcards):
    wanted_outputs = []
    samples_file = pd.read_csv(config['tsv_file'], sep="\t")
    comparaisons = get_all_pairwise_comparisons()
    token = config['token']
    wanted_outputs.append("results/" + config['token'] + "/ExecutionComplete.txt")
    for comparaison in comparaisons:
        if (config['filtering_activate'] == 'True'):
            wanted_outputs.append("results/{token}/in_house_method/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_in-house.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['gsea_activate'] == 'True'):
            wanted_outputs.append("results/{token}/GSEA-like/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_GSEA-like.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
            wanted_outputs.append("results/{token}/DESeq2/{treatment}_vs_{control}/{treatment}_vs_{control}_DESeq2_table.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['bagel_activate'] == 'True'):
            wanted_outputs.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.bf".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
            wanted_outputs.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_bagel_essentials_genes.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if config['mageck_mle_activate'] == 'True':
            wanted_outputs.append("results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if config['mageck_rra_activate'] == 'True':
            wanted_outputs.append("results/{token}/MAGeCK_RRA/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if config['crisphiermix_activate'] == 'True':
            wanted_outputs.append("results/{token}/CRISPhieRmix/{treatment}_vs_{control}/{treatment}_vs_{control}.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
    return wanted_outputs


rule MAGeCK_counts_normalize:
    input:
        config['count_table_file']
    output:
        config['normalized_count_table']
    params:
        name = "resources/" + config['token'] + "/screen",
    log:
        "logs/normalization/MAGeCK_counts_normalize.log"
    shell:
        "mageck count \
            -k {input} \
            -n {params.name} \
            --norm-method total > {log}"
