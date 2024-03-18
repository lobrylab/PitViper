import pandas as pd
import sys
import logging

# Get the input and output file names from the command line arguments
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Set up logging
log_file = snakemake.log[0]
logging.basicConfig(filename=log_file, level=logging.ERROR)

try:
    # Read the tab-separated file
    df = pd.read_csv(input_file, sep="\t", header=None)

    # Swap the columns
    df = df[[3, 0, 1, 2]]

    # Write the dataframe to a new file
    df.to_csv(output_file, sep="\t", header=False, index=False)

    # Add 5 empty lines at the beginning of the file
    with open(output_file, "r+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write("#\n" * 6 + content)

except Exception as e:
    logging.error("Caught an error: " + str(e))
