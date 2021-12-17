import click
import os
import subprocess


def run_pitviper(configfile, dry_run, jobs, notebook):
    if notebook != "":
        nb_opt = "--edit-notebook %s" % notebook
    else:
        nb_opt = ""
    if dry_run:
        cmd = "snakemake -s workflow/Snakefile -n --configfile {configfile} --use-conda --cores {jobs}".format(configfile=configfile, jobs=jobs)
    elif not dry_run:
        cmd = "snakemake -s workflow/Snakefile --configfile {configfile} --use-conda --cores {jobs} {nb_opt}".format(configfile=configfile, jobs=jobs, nb_opt=nb_opt)
    print('Command:', cmd)
    os.system(cmd)


@click.command()
@click.option('--configfile', help='Path to configuration file. Must be in YAML format.', type=str, required=True)
@click.option('--dry_run', help='If dry run is enabled, commands won\'t be executed. A message will be printed instead.', default=False, type=bool)
@click.option('--jobs', help='Define number of jobs to run.', default=1, type=int)
@click.option('--notebook', help='Output file(s) of a PitViper rule with "notebook" entry.', type=str, required=False, default="")
def main(configfile, dry_run, jobs, notebook):
    try:
        print('Working directory:', os.getcwd())
        print('Trying to open:', configfile)
        f = open(configfile, 'r')
        f.close()
    except IOError:
        print("Error: cannot open and read configuration file {configfile}".format(configfile=configfile))
    if dry_run:
        run_pitviper(configfile=configfile, dry_run=True, jobs=jobs, notebook=notebook)
    else:
        run_pitviper(configfile=configfile, dry_run=False, jobs=jobs, notebook=notebook)

if __name__ == '__main__':
    main()
