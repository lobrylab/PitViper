import click
import pandas as pd
import sys
import numpy as np

@click.command()
@click.option('--output', '-o', help='Path to write output file.', type=str, required=True)
@click.option('--count_table', '-c', help='Path to count table.', type=str, required=True)
@click.option('--samples', '-s', help='List of samples, separated by a comma.', type=str, required=True)
@click.option('--replicates', '-r', multiple=True)
def main(count_table, samples, replicates, output):
    # Open count table with pandas
    count_table = pd.read_table(count_table)

    samples_list = samples.split(',')

    if len(replicates) != len(samples_list):
        print(len(replicates), 'replicates lists and', len(samples_list), 'samples.')
        print('Error: can\'t match replicates with samples.')
        sys.exit(0)

    table = []

    for i in range(len(samples_list)):
        sample = samples_list[i]
        for j in range(len(replicates[i].split(','))):
            replicate = replicates[i].split(',')[j]
            # row = [ replicate, sample, j+1 ]
            row = [replicate, sample]
            table.append(row)

    table = pd.DataFrame(table, columns=['sample', 'replicate'])
    table.to_csv(output, index=False)

if __name__ == '__main__':
    main()
