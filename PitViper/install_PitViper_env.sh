#!/bin/bash -i

### Create a conda env from big_env.yaml file named 'pitviper_env'.
mamba env create -n pitviper_env -f pitviper_env.yaml

conda activate pitviper_env

### Install CRISPhieRmix R package in 'pitviper_env'.
Rscript -e "devtools::install_github('timydaley/CRISPhieRmix')"
