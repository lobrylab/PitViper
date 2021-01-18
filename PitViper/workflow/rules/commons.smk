def getTreatmentIds(wildcards):
    samples_file = pd.read_csv(config['inputs']['samples'], sep=",")
    list_ids = list(samples_file[wildcards.treatment].values)
    return ",".join(list_ids)


def getControlIds(wildcards):
    samples_file = pd.read_csv(config['inputs']['samples'], sep=",")
    list_ids =  list(samples_file[wildcards.control].values)
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

    samples_file = pd.read_csv(config['inputs']['samples'], sep=",")
    samples_list = list(samples_file.columns)
    comparaisons = pairwiseComparaisons(samples_list)


    l = [{'treatment': samples_list[duo[1]], 'control':samples_list[duo[0]]} for duo in comparaisons]

    return l

get_all_pairwise_comparaisons()




def get_pipeline_outputs(wildcards):
    wanted_outputs = []

    comparaisons = get_all_pairwise_comparaisons()

    print(comparaisons)

    token = config['token']

    for comparaison in comparaisons:
        wanted_outputs.append("results/{token}/BAGEL/{treatment}_vs_{control}/{treatment}_vs_{control}_BAGEL_output.bf".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        wanted_outputs.append("results/{token}/MAGeCK_MLE/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
        wanted_outputs.append("results/{token}/MAGeCK_RRA/{treatment}_vs_{control}/{treatment}_vs_{control}.gene_summary.txt".format(token = token, treatment = comparaison['treatment'], control = comparaison['control']))
    return wanted_outputs
