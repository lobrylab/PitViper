#!/bin/bash -i

### Activation of pitviper_env.
conda activate pitviper_env

### Run PitViper application with Flask.
(cd ./PitViper && python3 gui/waitress_server.py)
