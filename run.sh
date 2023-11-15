#!/bin/bash

# This script sets up the environment and runs the PitViper Flask application

# Exit if any command fails
set -e

# Constants for installer and mode
readonly CONDA_INSTALLER="conda"
readonly MAMBA_INSTALLER="mamba"
readonly YAML_MODE="yaml"
readonly LOCK_MODE="lock"

installer="mamba"
mode="yaml"
noflask="false"

# Parse the command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --help)
      echo "Usage: $0 [options]"
      echo "Options:"
      echo "  --mode           The mode to use for creating the conda environment (yaml or lock). Default is yaml."
      echo "  --installer      The installer to use for creating the conda environment (conda or mamba). Default is mamba."
      echo "  --noflask        Do not run the Flask application."
      echo "  --help           Display this help message."
      exit 0
      ;;
    --mode)
      mode="$2"
      shift 2
      ;;
    --installer)
      installer="$2"
      shift 2
      ;;
    --noflask)
      noflask="true"
      shift 1
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

# # Initialize installer variable before conditional checks
# installer="${1:-$MAMBA_INSTALLER}"
# mode="${2:-$YAML_MODE}"

# Set configurable variables using environment variables
app_env_var=${APP_ENV_VAR:-FLASK_APP}
conda_env_name=${CONDA_ENV_NAME:-pitviper_env}
app_path=${APP_PATH:-PitViper/gui/app.py}

# Check if installer is 'conda' or 'mamba'
if [[ "$installer" == "$CONDA_INSTALLER" ]]; then
    echo "Using conda installer"
elif [[ "$installer" == "$MAMBA_INSTALLER" ]]; then
    echo "Using mamba installer"
else
    # Set default installer to conda
    installer="$MAMBA_INSTALLER"
    echo "No installer specified. Using mamba installer by default."
fi

# Check if mode is 'yaml' or 'lock'
if [[ "$mode" == "$YAML_MODE" ]]; then
    echo "Using PitViper/environment.yaml file"
elif [[ "$mode" == "$LOCK_MODE" ]]; then
    echo "Using PitViper/environment_lock.yaml file"
else
    # Set default mode to yaml
    mode="$YAML_MODE"
    echo "No mode specified. Using PitViper/environment.yaml file by default."
fi


# Check that the required file exists
test -e "$app_path" || { echo "Error: $app_path does not exist"; exit 1; }

# Use absolute path for file
app_full_path=$(realpath "$app_path")

# Create a conda environment for the application
if conda env list | grep -qsw "$conda_env_name"; then
    echo "Activating conda environment: $conda_env_name"
    conda_path=$(conda info | grep -i 'base environment' | awk '{print $4}')

    # Activate conda environment
    source $conda_path/etc/profile.d/conda.sh
    conda activate "$conda_env_name"
else

    # Create conda environment
    echo "Creating conda environment: $conda_env_name"

    # Check if installer is 'conda' or 'mamba'
    if [[ "$installer" == "$CONDA_INSTALLER" ]]; then

        # Install packages using conda
        echo "Installing packages using conda..."
        if [[ "$mode" == "$YAML_MODE" ]]; then

            # Install packages using conda and environment.yaml file
            conda env create -f PitViper/environment.yaml -n "$conda_env_name"
        elif [[ "$mode" == "$LOCK_MODE" ]]; then

            # Install packages using conda and environment.lock file
            conda create --file PitViper/environment_lock.yaml -n "$conda_env_name"
        fi
    elif [[ "$installer" == "$MAMBA_INSTALLER" ]]; then

        # Install packages using mamba
        echo "Installing packages using mamba..."

        # Check if mamba package is already installed
        if conda list -n base | grep -qws "mamba"; then
            echo "Mamba package is already installed."

        # Install mamba package if it's not already installed
        else
            echo "Installing Mamba package..."
            conda install -y -c conda-forge mamba
        fi

        # Check if mode is 'yaml' or 'lock'
        if [[ "$mode" == "$YAML_MODE" ]]; then

            # Install packages using mamba and environment.yaml file
            mamba env create -f PitViper/environment.yaml -n "$conda_env_name"
        elif [[ "$mode" == "lock" ]]; then

            # Install packages using mamba and environment.lock file
            mamba create --file PitViper/environment_lock.yaml -n "$conda_env_name"
        fi
    fi
    # Activate conda environment
    conda_path=$(conda info | grep -i 'base environment' | awk '{print $4}')
    source $conda_path/etc/profile.d/conda.sh
    conda activate "$conda_env_name"

    # Install CRISPhieRmix
    curl -L https://github.com/lobrylab/CRISPhieRmix/tarball/e400f21 -o PitViper/CRISPhieRmix.tar.gz
    Rscript -e 'install.packages("PitViper/CRISPhieRmix.tar.gz", repos = NULL, type="source")'
    rm -rf PitViper/CRISPhieRmix.tar.gz
    
    # Install BAGEL2
    # Check if PitViper/workflow/scripts/bagel already exists, if not clone the repository
    bagel_dir="PitViper/workflow/scripts/bagel"
    if [[ ! -d "$bagel_dir" ]]; then
        git clone https://github.com/lobrylab/bagel.git "$bagel_dir"
    fi

    # Freeze the conda environment to a YAML file
    conda env export > PitViper/environment_freeze.yaml
fi

# Set the Flask app environment variable if it's not already set
if [[ "${!app_env_var}" != "$app_full_path" ]]; then
    export "$app_env_var=$app_full_path"
fi

# If noflask is set to false, run the Flask application
if [[ "$noflask" == "false" ]]; then
    # Run the Flask application
    (cd PitViper/ && flask run)
fi
