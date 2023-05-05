"""
PitViper CLI: a pipeline for the analysis of single-cell RNA-seq data

Usage:
    pitviper.py --configfile <configfile> [--dry_run] [--jobs <jobs>] [--notebook <notebook>]
    
Options:
    --configfile <configfile>   Path to configuration file. Must be in YAML format.
    --dry_run                   If dry run is enabled, commands won't be executed. A message will be printed instead.
    --jobs <jobs>               Define number of Snakemake rules to run in parallel. [default: 1]
    --notebook <notebook>       Output file(s) of a PitViper rule with "notebook" entry.
    -h, --help                  Show this screen.
    --version                   Show version.
    
Examples:
    pitviper.py --configfile config.yaml
    pitviper.py --configfile config.yaml --dry_run
    pitviper.py --configfile config.yaml --jobs 4
    pitviper.py --configfile config.yaml --notebook notebook.html
"""
import os

import click
import yaml


def _check_notebook_arguments(notebook):
    # Get notebook location
    notebook_location = os.path.dirname(notebook)
    # Check if notebook location exists
    if not os.path.exists(notebook_location):
        raise IOError(f"Error: notebook location `{notebook_location}` does not exist.")
    # Check if notebook location is a directory
    if not os.path.isdir(notebook_location):
        raise IOError(
            f"Error: notebook location {notebook_location} is not a directory."
        )
    # Check if notebook location has write permissions
    if not os.access(notebook_location, os.W_OK):
        raise IOError(
            f"Error: notebook location {notebook_location} does not have write permissions."
        )
    # Check if notebook already exists
    if os.path.exists(notebook):
        raise IOError(f"Error: notebook {notebook} already exists.")
    # Check if notebook extension is .ipynb
    if not notebook.endswith(".ipynb"):
        raise IOError(f"Error: notebook {notebook} does not have .ipynb extension.")


def _check_jobs_argument(jobs):
    try:
        jobs = int(jobs)
    except ValueError as exc:
        raise IOError("Error: jobs argument must be an integer.") from exc
    if jobs < 1:
        raise IOError("Error: jobs argument must be greater than 0.")
    return jobs


def run_pitviper(configfile, dry_run, jobs, notebook):
    """
    Run PitViper pipeline.
    """
    if notebook != "":
        _check_notebook_arguments(notebook)
        nb_opt = f"--edit-notebook {notebook}"
    else:
        nb_opt = ""
    if dry_run:
        cmd = f"snakemake -s workflow/Snakefile -n --configfile {configfile} --use-conda --cores {jobs}"
    elif not dry_run:
        cmd = f"snakemake -s workflow/Snakefile --configfile {configfile} --use-conda --cores {jobs} {nb_opt}"
    print("Command:", cmd)
    os.system(cmd)


def _check_file_entries(file):
    """
    Check if configuration file contains required entries.
    """
    file = open(file, "r", encoding="utf-8")
    config = yaml.safe_load(file)
    missing_entries = []
    required_entries = [
        "bagel_activate",
        "bed_annotation_file",
        "bowtie_activate",
        "bowtie_mapping_method",
        "controls_file",
        "count_table_file",
        "counts_threshold",
        "crisphiermix_activate",
        "directional_scoring_method_activate",
        "directional_scoring_method_fdr_threshold",
        "directional_scoring_method_guides_threshold",
        "directional_scoring_method_log2_threshold",
        "essentials",
        "genome_version_rose",
        "jobs",
        "length_3_adapter",
        "length_5_adapter",
        "library_file",
        "mageck_count_N",
        "mageck_count_activate",
        "mageck_count_all_align",
        "mageck_count_length",
        "mageck_count_normalization",
        "mageck_count_rev_comp",
        "mageck_mle_activate",
        "mageck_mle_adj",
        "mageck_mle_mean_var",
        "mageck_mle_normalization",
        "mageck_mle_outliers",
        "mageck_mle_perm_N",
        "mageck_mle_perm_all",
        "mageck_rra_LFC",
        "mageck_rra_activate",
        "mageck_rra_adj",
        "mageck_rra_count_min",
        "mageck_rra_criteria",
        "mageck_rra_normalization",
        "mageck_rra_pthreshold",
        "mageck_rra_remove",
        "nonessentials",
        "normalized_count_table",
        "screen_type",
        "ssrea_activate",
        "ssrea_ranking_method",
        "start_from",
        "token",
        "tsv_file",
    ]
    for entry in required_entries:
        if entry not in config:
            missing_entries.append(entry)
    if len(missing_entries) > 0:
        raise IOError(
            f"Error: configuration file is missing required entries: {missing_entries}"
        )


def _check_file_format(file):
    """
    Check if file is in YAML format.
    """
    try:
        yaml.safe_load(file)
    except yaml.YAMLError as exc:
        raise IOError(f"Error: file {file} is not in YAML format.") from exc


def _check_file_permissions(file):
    """
    Check if file has read permissions.
    """
    if not os.access(file, os.R_OK):
        raise IOError(f"Error: file {file} does not have read permissions.")


def _check_file_size(file):
    """
    Check if file size is greater than 0.
    """
    if os.stat(file).st_size == 0:
        raise IOError(f"Error: file {file} is empty.")


def _check_file_type(file):
    """
    Check if file is a file.
    """
    if not os.path.isfile(file):
        raise IOError(f"Error: file {file} is not a file or does not exist.")


def _check_configuration(file: str):
    """
    Check if file exists and has read permissions.
    """
    _check_file_type(file)
    _check_file_permissions(file)
    _check_file_size(file)
    _check_file_format(file)
    _check_file_entries(file)
    print("Configuration file is valid.")


@click.command()
@click.option(
    "--configfile",
    help="Path to configuration file. Must be in YAML format.",
    type=str,
    required=True,
)
@click.option(
    "--dry_run",
    help="If dry run is enabled, commands won't be executed. A message will be printed instead.",
    default=False,
    type=bool,
)
@click.option(
    "--jobs",
    help="Define number of Snakemake rules to run in parallel.",
    default=1,
    type=int,
)
@click.option(
    "--notebook",
    help='Output file(s) of a PitViper rule with "notebook" entry.',
    type=str,
    required=False,
    default="",
)
def main(configfile, dry_run, jobs, notebook):
    """
    PitViper CLI: a pipeline for the analysis of single-cell RNA-seq data.
    """
    # Check if configuration file is valid
    _check_configuration(configfile)
    jobs = _check_jobs_argument(jobs)
    run_pitviper(configfile=configfile, dry_run=dry_run, jobs=jobs, notebook=notebook)


if __name__ == "__main__":
    main()
