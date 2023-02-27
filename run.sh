#!/bin/bash

# This script sets up the environment and runs the PitViper Flask application

# Exit if any command fails
set -e

# Set configurable variables using environment variables
app_env_var=${APP_ENV_VAR:-FLASK_APP}
conda_env_name=${CONDA_ENV_NAME:-pitviper_env}
app_path=${APP_PATH:-PitViper/gui/app.py}

# Check that the required file
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
    echo "Creating conda environment: $conda_env_name"
    conda install -c conda-forge -y mamba
    mamba env create -f PitViper/environment.yaml -n "$conda_env_name"
    conda_path=$(conda info | grep -i 'base environment' | awk '{print $4}')
    source $conda_path/etc/profile.d/conda.sh
    conda activate "$conda_env_name"
    Rscript -e 'install.packages("PitViper/workflow/scripts/CRISPhieRmix-1.1.tar.gz", repos = NULL, type="source")'
fi

# Set the Flask app environment variable if it's not already set
if [[ "${!app_env_var}" != "$app_full_path" ]]; then
    export "$app_env_var=$app_full_path"
fi

(cd PitViper/ && flask run)
