#!/bin/bash

### Activation of pitviper_env.
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh
conda activate pitviper_env


### Get absolute path of Flask application.
app_full_path=$(realpath PitViper/gui/app.py)


### Change $FLASK_APP value if not equal to app absolute path.
if [ "$FLASK_APP" != "$app_full_path" ]
then
  export FLASK_APP="$app_full_path"
fi


### Run PitViper application with Flask.
(cd PitViper/ && flask run)


# Open the jupyter Notebook web-page in results/ directory.
jupyter notebook PitViper/results/

