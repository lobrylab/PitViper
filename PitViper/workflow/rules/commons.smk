
def getFastqFiles(wildcards):
    """Return a list of experiment's files  in fastq format."""
    samples_sheet = pd.read_csv(config['tsv_file'], sep="\t")
    if config['mageck_count_activate'] == 'True':
        fastqs = samples_sheet.fastq.values
    else:
        fastqs = []
    return fastqs


def getBamFiles(wildcards):
    """Return a list of experiment's files in bam format."""
    samples_sheet = pd.read_csv(config['tsv_file'], sep="\t")
    if config['mageck_count_activate'] == 'True':
        bams = samples_sheet.bam.values
    elif config['bowtie_activate'] == 'True':
         bams = []
         replicates = samples_sheet.replicate.values
         for replicate in replicates:
             bam = "results/%s/bam/%s.bam" % (config['token'], replicate)
             bams.append(bam)
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
    if config['mageck_count_activate'] == 'True':
        if 'fastq' in samples_sheet.columns:
            fastqs = []
            for i in range(0, len(samples_sheet.index)):
                fastqs.append(samples_sheet.loc[i].fastq)
            return " ".join(fastqs)
        elif 'bam' in samples_sheet.columns:
            bam_files = []
            for i in range(0, len(samples_sheet.index)):
                bam_files.append(samples_sheet.loc[i].bam)
            return " ".join(bam_files)
    elif config['bowtie_activate'] == 'True':
        bams = []
        replicates = samples_sheet.replicate.values
        for replicate in replicates:
            bam = "results/%s/bam/%s.bam" % (config['token'], replicate)
            bams.append(bam)
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


def bagel_bf_columns(wildcards):
    """Return a list of column numbers for BAGEL output file. Use in bagel.smk"""
    file = "results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL.foldchange".format(token=wildcards.token, treatment=wildcards.treatment, control=wildcards.control)
    content = pd.read_csv(file, sep="\t")
    out = ",".join([ str(i + 1) for i in range(len(content.columns)-2)])
    return out


def generatedResults(wildcards):
    """Return a list of all the results files to integrate in the report. Use in integration.smk"""
    results = []
    comparaisons = get_all_pairwise_comparisons()
    token = config['token']
    for comparaison in comparaisons:
        if (config['bagel_activate'] == 'True'):
            results.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.bf".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
            results.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.pr".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
            results.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL1_output.bf".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
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


def get_pipeline_outputs(wildcards):
    wanted_outputs = []
    samples_file = pd.read_csv(config['tsv_file'], sep="\t")
    comparaisons = get_all_pairwise_comparisons()
    token = config['token']
    wanted_outputs.append(f"results/{config['token']}/Report.ipynb")
    # If `fastq` or `bam` columns are present in the tsv file, we add "results/{token}/fastqc" to the list of outputs
    if 'fastq' in samples_file.columns or 'bam' in samples_file.columns:
        wanted_outputs.append(f"results/{config['token']}/fastqc")
    for comparaison in comparaisons:
        if (config['directional_scoring_method_activate'] == 'True'):
            wanted_outputs.append("results/{token}/directional_scoring_method/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_directional_scoring_method.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['ssrea_activate'] == 'True'):
            wanted_outputs.append("results/{token}/SSREA/{treatment}_vs_{control}/{treatment}_vs_{control}_all-elements_SSREA.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
            wanted_outputs.append("results/{token}/DESeq2/{treatment}_vs_{control}/{treatment}_vs_{control}_DESeq2_table.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if (config['bagel_activate'] == 'True'):
            wanted_outputs.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.bf".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
            wanted_outputs.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.pr".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if config['mageck_mle_activate'] == 'True':
            wanted_outputs.append("results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if config['mageck_rra_activate'] == 'True':
            wanted_outputs.append("results/{token}/MAGeCK_RRA/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if config['crisphiermix_activate'] == 'True':
            wanted_outputs.append("results/{token}/CRISPhieRmix/{treatment}_vs_{control}/{treatment}_vs_{control}.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if config['bed_annotation_file'] != "":
            wanted_outputs.append("resources/{token}/annotation_ROSE_REGION_TO_GENE.txt".format(token = token))
    return wanted_outputs
