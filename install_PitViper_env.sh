#!/bin/bash

# This script sets up the environment and runs the PitViper Flask application

# Get first argument
installer=$1

# Check if installer is 'conda' or 'mamba'
if [[ "$installer" == "conda" ]]; then
    echo "Using conda installer"
elif [[ "$installer" == "mamba" ]]; then
    echo "Using mamba installer"
else
    # Set default installer to conda
    installer="mamba"
    echo "No installer specified. Using mamba installer by default."
fi

# Exit if any command fails
set -e

# Set configurable variables using environment variables
app_env_var=${APP_ENV_VAR:-FLASK_APP}
conda_env_name=${CONDA_ENV_NAME:-pitviper_env}
app_path=${APP_PATH:-PitViper/gui/app.py}

# Check that the required file exists
test -e "$app_path" || { echo "Error: $app_path does not exist"; exit 1; }

# Use absolute path for file
app_full_path=$(realpath "$app_path")

# Create a conda environment for the application
if conda env list | grep -qsw "$conda_env_name"; then
    echo "Activating conda environment: $conda_env_name"
    conda_path=$(conda info | grep -i 'base environment' | awk '{print $4}')
    source $conda_path/etc/profile.d/conda.sh
    conda activate "$conda_env_name"
else
    # Create conda environment
    echo "Creating conda environment: $conda_env_name"
    # Check if installer is 'conda' or 'mamba'
    if [[ "$installer" == "conda" ]]; then
        # Install packages using conda
        echo "Installing packages using conda..."
        conda env create -f PitViper/environment.yaml -n "$conda_env_name"
    elif [[ "$installer" == "mamba" ]]; then
        # Install packages using mamba
        echo "Installing packages using mamba..."

        # Check if mamba package is already installed
        if conda list -n base | grep -q -w "mamba"; then
            echo "Mamba package is already installed."
        # Install mamba package if it's not already installed
        else
            echo "Installing Mamba package..."
            conda remove -y libarchive && conda install -y -c conda-forge libarchive=3.6.2 mamba
        fi
        mamba env create -f PitViper/environment.yaml -n "$conda_env_name"
    fi
    conda_path=$(conda info | grep -i 'base environment' | awk '{print $4}')
    source $conda_path/etc/profile.d/conda.sh
    conda activate "$conda_env_name"
    curl -L https://github.com/timydaley/CRISPhieRmix/tarball/master -o /PitViper/CRISPhieRmix.tar.gz
    #Rscript -e 'install.packages(c("Rcpp", "statmod", "sn", "nloptr"), repos = "CRAN")'
    Rscript -e 'install.packages("/PitViper/CRISPhieRmix.tar.gz", repos = NULL, type="source")'
    rm -rf /PitViper/CRISPhieRmix.tar.gz
    #Rscript -e 'devtools::install_github("timydaley/CRISPhieRmix")'
    git clone https://github.com/hart-lab/bagel.git PitViper/workflow/scripts/bagel

    # Freeze the conda environment to a YAML file
    conda env export > PitViper/environment_freeze.yaml
fi