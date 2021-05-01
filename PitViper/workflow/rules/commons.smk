def getFastqFiles(wildcards):
    """Return a list of experiment's files  in fastq format."""
    samples_sheet = pd.read_csv(config['inputs']['tsv'], sep="\t")
    fastqs = samples_sheet.fastq.values
    return fastqs


def getBamFiles(wildcards):
    """Return a list of experiment's files in bam format."""
    samples_sheet = pd.read_csv(config['inputs']['tsv'], sep="\t")
    bams = samples_sheet.bam.values
    return bams


def getLabels(wildcards):
    """Return concatenation of experiment's labels. Needed for MAGeCK count software."""
    samples_sheet = pd.read_csv(config['inputs']['tsv'], sep="\t")
    labels = samples_sheet.replicate.values
    labels_str = ",".join(labels)
    return labels_str


def getFiles(wildcards):
    """Return concatenation of experiment's files. Needed for MAGeCK count software."""
    samples_sheet = pd.read_csv(config['inputs']['tsv'], sep="\t")
    if 'fastq' in samples_sheet.columns:
        fastqs = samples_sheet.fastq.values
        fastqs_str = " ".join(fastqs)
        return fastqs_str
    elif 'bam' in samples_sheet.columns:
        bams = samples_sheet.bam.values
        bams_str = " ".join(bams)
        return bams_str



def getTreatmentIdsLen(wildcards):
    samples_file = pd.read_csv(config['inputs']['tsv'], sep="\t")
    list_ids = samples_file.loc[samples_file["condition"] == wildcards.treatment].values
    return len(list_ids)


def getControlIdsLen(wildcards):
    samples_file = pd.read_csv(config['inputs']['tsv'], sep="\t")
    list_ids = samples_file.loc[samples_file["condition"] == wildcards.control].values
    return len(list_ids)


def getTreatmentIds(wildcards):
    samples_file = pd.read_csv(config['inputs']['tsv'], sep="\t")
    list_ids = samples_file.loc[samples_file.condition == wildcards.treatment].replicate.values
    return ",".join(list_ids)


def getControlIds(wildcards):
    samples_file = pd.read_csv(config['inputs']['tsv'], sep="\t")
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

    samples_file = pd.read_csv(config['inputs']['tsv'], sep="\t")
    samples_list = list(set(samples_file.condition.values))
    comparaisons = pairwiseComparaisons(samples_list)


    l = [{'treatment': samples_list[duo[1]], 'control':samples_list[duo[0]]} for duo in comparaisons]

    return l

def get_pipeline_outputs(wildcards):
    wanted_outputs = []

    comparaisons = get_all_pairwise_comparaisons()

    token = config['token']

    #wanted_outputs.append("results/" + config['token'] + "/reports/PitViper_report.ipynb")

    for comparaison in comparaisons:
        if (config['BAGEL']['activate'] == True):
            wanted_outputs.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.bf".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
            wanted_outputs.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_bagel_essentials_genes.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if config['MAGeCK']['MLE']['activate'] == 'True':
            wanted_outputs.append("results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if config['MAGeCK']['RRA']['activate'] == 'True':
            wanted_outputs.append("results/{token}/MAGeCK_RRA/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        if config['CRISPhieRmix']['activate'] == 'True':
            wanted_outputs.append("results/{token}/CRISPhieRmix/{treatment}_vs_{control}/{treatment}_vs_{control}.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))

    return wanted_outputs
