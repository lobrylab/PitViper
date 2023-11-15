import pandas as pd

# Handle Snakemake arguments.
cts_file = snakemake.input[0]
res_file = snakemake.input[1]
treatment = snakemake.params[0]
baseline = snakemake.params[1]
method = snakemake.params[2]
design_file = snakemake.input[2]

# Process counts file (must be tab delimited and have a 'sgRNA' column).
cts = pd.read_csv(cts_file, sep="\t", index_col="sgRNA")

# Read design file for replicate/condition associations.
design = pd.read_csv(design_file, sep="\t")

# Get treatment and baseline replicates.
treatment_reps = design[design["condition"] == treatment]["replicate"].tolist()
baseline_reps = design[design["condition"] == baseline]["replicate"].tolist()

# Filter counts to treatment and baseline replicates in one cts dataframe.
cts = cts[treatment_reps + baseline_reps]

# Count the number of replicates for each condition.
n_treatment = len(treatment_reps)
n_baseline = len(baseline_reps)

if method == "signal_to_noise" and n_treatment > 2 and n_baseline > 2:
    print("Using signal to noise method.")
    # Compute mean and standard deviation for each sgRNA for each condition.
    treatment_mean = cts[treatment_reps].mean(axis=1)
    treatment_std = cts[treatment_reps].std(axis=1)
    baseline_mean = cts[baseline_reps].mean(axis=1)
    baseline_std = cts[baseline_reps].std(axis=1)

    # Compute a signal to noise ratio for each sgRNA.
    table = (treatment_mean - baseline_mean) / (treatment_std + baseline_std)
    table.sort_values(ascending=False, inplace=True)

    # Add index as a column.
    table = table.reset_index()

    # Rename columns.
    table.columns = ["sgRNA", "ranking"]

    # Filter out sgRNAs with NaN values.
    table = table.dropna()

    # Reorder columns to signal_to_noise, sgRNA.
    table = table[["ranking", "sgRNA"]]
else:
    print("Using the Log2 of classes method.")
    # Open results file.
    table = pd.read_csv(res_file, sep="\t")

    # Rename columns 'log2FoldChange' to 'ranking'
    table.rename(columns={"log2FoldChange": "ranking"}, inplace=True)
    table.sort_values(by="ranking", ascending=False, inplace=True)

    # Filter out sgRNAs with NaN values.
    table = table.dropna()

    # Select only the 'ranking' and 'sgRNA' columns.
    table = table[["ranking", "sgRNA"]]

# Write results to file.
table.to_csv(snakemake.output[0], sep="\t", index=False)
