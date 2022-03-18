#!/bin/bash

### Activation of pitviper_env.
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh
conda activate pitviper_env

### Get absolute path of Flask application.
app_full_path=$(realpath gui/app.py)

### Show app.py absolute path.
echo $app_full_path

### Change $FLASK_APP value if not equal to app absolute path.
if [ "$FLASK_APP" != "$app_full_path" ]
then
  export FLASK_APP="$app_full_path"
fi

### Show new $FLASK_APP value.
echo $FLASK_APP

### Run PitViper application with Flask.
flask run --host=0.0.0.0

jupyter notebook --allow-root --no-browser --ip=0.0.0.0 results/
