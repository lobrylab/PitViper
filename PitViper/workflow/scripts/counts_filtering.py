import pandas as pd


def get_all_pairwise_comparisons(design):
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


def read_cts(cts_path):
    return pd.read_csv(cts_path, sep='\t', header=0)


def read_design(design_path):
    return pd.read_csv(design_path, sep="\t")

design = read_design(snakemake.params[0])
pairwise_comparisons = get_all_pairwise_comparisons(design)
cts = read_cts(snakemake.input[0])
cts = pd.melt(cts, id_vars=['sgRNA', 'Gene'], value_vars=cts.columns.values[2:])

cts = pd.merge(cts, design, left_on='variable', right_on='replicate')[["sgRNA", "Gene", "value", "condition", "replicate"]]

cts = cts.groupby(['sgRNA', 'Gene', 'condition']).median().reset_index()


for comparison in pairwise_comparisons:
    label = comparison['treatment'] + "_vs_" + comparison['control']
    cts['below_threshold'] = cts['value'] < 100


matrix = pd.pivot_table(cts, values='below_threshold', index=['sgRNA', 'Gene'], columns=['condition']).reset_index()
matrix['keep'] = ~matrix.all(axis='columns', bool_only = bool)

guides_to_keep = matrix.loc[matrix['keep'] == True].sgRNA.values

print("Len:", len(guides_to_keep))

cts_normlized = read_cts(snakemake.input[0])
cts_normlized = cts_normlized.loc[cts_normlized.sgRNA.isin(guides_to_keep)]
cts_normlized.to_csv(snakemake.output[0], sep = "\t", index = False)

cts_raw = read_cts(snakemake.input[1])
cts_raw = cts_raw.loc[cts_raw.sgRNA.isin(guides_to_keep)]
cts_raw.to_csv(snakemake.output[1], sep = "\t", index = False)
