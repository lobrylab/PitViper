import click
import os

def run_pitviper(configfile, dry_run):
    if dry_run:
        cmd = "snakemake -d PitViper/ -s PitViper/workflow/Snakefile -n --configfile {configfile} --use-conda --cores 1".format(configfile=configfile)
    elif not dry_run:
        cmd = "snakemake -d PitViper/ -s PitViper/workflow/Snakefile --configfile {configfile} --use-conda --cores 1".format(configfile=configfile)
    print('Command:', cmd)
    os.system(cmd)


@click.command()
@click.option('--run_snakemake', help='If True, snakemake pipeline will be executed with choosen parameters.', default=False, type=bool, required=True)
@click.option('--configfile', help='Path to configuration file. Must be in YAML format.', type=str, required=True)
@click.option('--dry_run', help='If dry run is enabled, commands won\'t be executed. A message will be printed instead.', default=False, type=bool)
def main(run_snakemake, configfile, dry_run):
    try:
        configfile = 'PitViper/{configfile}'.format(configfile=configfile)
        f = open(configfile)
        if dry_run:
            run_pitviper(configfile=configfile, dry_run=True)
            f.close()
        elif run_snakemake:
            run_pitviper(configfile=configfile, dry_run=False)
            f.close()
        else:
            print('Dry run activated or snakemake not activated, pipeline won\'t be executed.')
    except IOError:
        print("Error: cannot open and read configuration file {configfile}".format(configfile=configfile))


if __name__ == '__main__':
    main()
