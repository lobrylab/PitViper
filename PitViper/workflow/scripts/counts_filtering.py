import pandas as pd

def main(threshold=0):
    """
    Main counts filtering function.
    Filter sgRNA with all conditions having a median count value below threshold.
    """
    # Read design table.
    design = pd.read_csv(snakemake.params[0], sep="\t")
    # Read count table.
    cts = pd.read_csv(snakemake.input[0], sep='\t', header=0)
    # Unpivot the dataframe from wide to long format.
    cts = pd.melt(cts, id_vars=['sgRNA', 'Gene'],
                  value_vars=cts.columns.values[2:])
    # Merge design and count table
    cts = pd.merge(cts, design,
                   left_on='variable',
                   right_on='replicate')[["sgRNA",
                                          "Gene",
                                          "value",
                                          "condition",
                                          "replicate"]]
    # Compute median count value by sgRNA for each condition.
    cts = cts.groupby(['sgRNA', 'Gene', 'condition']).median().reset_index()
    # Mark sgRNA with median counts value below threshold.
    cts['to_remove'] = cts['value'] < threshold
    # Pivot dataframe to create a boolean dataframe
    cts = pd.pivot_table(cts, values='to_remove',
                         index=['sgRNA', 'Gene'],
                         columns=['condition']).reset_index()
    # Filter sgRNA with all conditions having a median count value below threshold.
    cts = cts[~cts.all(axis='columns', bool_only=bool)]
    # Get name of guides to keep.
    guides_to_keep = cts.sgRNA.values

    cts_normlized = pd.read_csv(snakemake.input[0], sep='\t', header=0)
    cts_normlized = cts_normlized.loc[cts_normlized.sgRNA.isin(guides_to_keep)]
    cts_normlized.to_csv(snakemake.output[0], sep="\t", index=False)

    cts_raw = pd.read_csv(snakemake.input[1], sep='\t', header=0)
    cts_raw = cts_raw.loc[cts_raw.sgRNA.isin(guides_to_keep)]
    cts_raw.to_csv(snakemake.output[1], sep="\t", index=False)    


if __name__ == '__main__':
    main(threshold=int(snakemake.params[1]))
