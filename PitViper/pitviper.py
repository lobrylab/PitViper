import click
import os
import subprocess


def run_pitviper(configfile, dry_run, jobs):
    if dry_run:
        cmd = "snakemake -s workflow/Snakefile -n --configfile {configfile} --use-conda --cores 1".format(configfile=configfile)
    elif not dry_run:
        cmd = "snakemake -s workflow/Snakefile --configfile {configfile} --use-conda --cores {jobs}".format(configfile=configfile, jobs=jobs)
    print('Command:', cmd)
    os.system(cmd)


@click.command()
@click.option('--run_snakemake', help='If True, snakemake pipeline will be executed with choosen parameters.', default=False, type=bool, required=True)
@click.option('--configfile', help='Path to configuration file. Must be in YAML format.', type=str, required=True)
@click.option('--dry_run', help='If dry run is enabled, commands won\'t be executed. A message will be printed instead.', default=False, type=bool)
@click.option('--jobs', help='Define number of jobs to run.', default=1, type=int)
def main(run_snakemake, configfile, dry_run, jobs):
    try:
        print('Working directory:', os.getcwd())
        print('Trying to open:', configfile)
        f = open(configfile, 'r')
        f.close()
    except IOError:
        print("Error: cannot open and read configuration file {configfile}".format(configfile=configfile))
    if dry_run:
        run_pitviper(configfile=configfile, dry_run=True)
    elif run_snakemake:
        run_pitviper(configfile=configfile, dry_run=False, jobs=jobs)
    else:
        print('Dry run activated or snakemake not activated, pipeline won\'t be executed.')


if __name__ == '__main__':
    main()
