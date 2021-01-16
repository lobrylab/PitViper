



def get_pipeline_outputs(wildcards):
    wanted_outputs = []

    samples = [str(i) for i in range(5)]

    samples = ["test"]

    #wanted_outputs.extend(list(map(lambda name: f"results/development/MAGeCK_MLE/{name}.gene_summary.txt", samples)))

#    wanted_outputs.extend(list(map(lambda name: f"results/development/{name}_succeed.txt", samples)))
    wanted_outputs.extend(list(map(lambda name: f"results/development/notebooks/MAGeCK_MLE_{name}.txt", samples)))

    return wanted_outputs
