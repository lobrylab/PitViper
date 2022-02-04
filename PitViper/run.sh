#!/bin/bash -i

### Activation of pitviper_env.
conda activate pitviper_test

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
flask run

jupyter notebook results/
