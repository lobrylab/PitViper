#!/bin/bash

### Activation of pitviper_env.
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh
conda activate pitviper_env

jupyter notebook --allow-root --no-browser --ip=0.0.0.0 PitViper/PitViper/results/

