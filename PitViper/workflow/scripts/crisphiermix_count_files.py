import pandas as pd
from functools import reduce
import click


def get_control(table, samples_table, control):
    header = ["sgRNA", "Gene"]
    replicates = samples_table.loc[samples_table.condition == control][
        "replicate"
    ].values
    header.extend(replicates)
    replicates_counts = table[header]
    return replicates_counts


def get_treatment(table, samples_table, treatment):
    header = ["sgRNA", "Gene"]
    replicates = samples_table.loc[samples_table.condition == treatment][
        "replicate"
    ].values
    header.extend(replicates)
    replicates_counts = table[header]
    return replicates_counts


@click.command()
@click.option("--file", default=1, help="Samples File.", required=True, type=str)
@click.option("--counts", default=1, help="Counts File.", required=True, type=str)
@click.option(
    "--directory", default=1, help="Output directory.", required=True, type=str
)
@click.option("--control", default=1, help="Control name.", required=True, type=str)
@click.option("--treatment", default=1, help="Treatment name.", required=True, type=str)
@click.option("--dryrun", default=False, help="Dry run.", type=bool)
def main(file, counts, control, treatment, directory, dryrun):
    """
    Merge the counts from the control and treatment samples into a single file.
    """
    table = pd.read_csv(counts, sep="\t")
    samples_table = pd.read_csv(file, sep=",")

    cont = get_control(table, samples_table, control)
    trea = get_treatment(table, samples_table, treatment)

    data_frames = [cont, trea]

    df_merged = reduce(
        lambda left, right: pd.merge(left, right, on=["sgRNA", "Gene"], how="outer"),
        data_frames,
    )

    file_name = directory + f"{treatment}_vs_{control}_count_matrix.txt"

    if not dryrun:
        df_merged.to_csv(file_name, index=False, sep="\t")
    else:
        print(file_name)
        print(df_merged)


if __name__ == "__main__":
    main()
