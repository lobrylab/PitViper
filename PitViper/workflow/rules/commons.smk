from natsort import natsorted

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
    labels = list(set(samples_sheet.replicate.values))
    labels_str = ",".join(labels)
    return labels_str


def getFiles(wildcards):
    """Return concatenation of experiment's files. Needed for MAGeCK count software."""
    samples_sheet = pd.read_csv(config['tsv_file'], sep="\t")
    if 'fastq' in samples_sheet.columns:
        conditions = list(set(samples_sheet.condition.values))
        fastqs = ""
        for condition in conditions:
            fastqs_condition = samples_sheet.loc[samples_sheet.condition == condition].fastq.values
            fastqs = fastqs + " " + " ".join(fastqs_condition)
        return fastqs
    elif 'bam' in samples_sheet.columns:
        conditions = list(set(samples_sheet.condition.values))
        bams = ""
        for condition in conditions:
            bams_condition = samples_sheet.loc[samples_sheet.condition == condition].bam.values
            bams = bams + " " + " ".join(bams_condition)
        return bams



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


def get_all_pairwise_comparaisons():
    def compareTuples(t1, t2):
        comparaison = set(t1) & set(t2)
        return (len(comparaison) == len(t1)) & (len(comparaison) == len(t2))


    def tupleInList(t, comparaisons):
        inList = False
        for comparaison in comparaisons:
            if compareTuples(t, comparaison):
                inList = True
                break
        return inList


    def pairwiseComparaisons(samples_list):
        comparaisons = []
        for i in range(0, len(samples_list)):
            for j in range(0, len(samples_list)):
                if (i != j) and not (tupleInList((i, j), comparaisons)):
                    comparaisons.append((i, j))
        return comparaisons

    samples_file = pd.read_csv(config['tsv_file'], sep="\t")
    samples_list = natsorted(list(set(samples_file.condition.values)))
    comparaisons = pairwiseComparaisons(samples_list)

    l = [{'treatment': samples_list[duo[1]], 'control':samples_list[duo[0]]} for duo in comparaisons]
    return l

def get_pipeline_outputs(wildcards):
    wanted_outputs = []

    comparaisons = get_all_pairwise_comparaisons()

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
