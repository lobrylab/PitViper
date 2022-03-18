#!/bin/bash

### Create a conda env from pitviper_env.yaml file named 'pitviper_env'.
mamba env create -n pitviper_env -f PitViper/environment.yaml

### Activation of pitviper_env.
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh
conda activate pitviper_env


### Install CRISPhieRmix R package in 'pitviper_env'.
Rscript -e 'install.packages("PitViper/workflow/scripts/CRISPhieRmix-1.1.tar.gz", repos = NULL, type="source")'
