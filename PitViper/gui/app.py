#!/home/paularthur/miniconda3/bin/python
import os
from flask import Flask, render_template, request, redirect, url_for
from werkzeug.utils import secure_filename

UPLOAD_FOLDER = 'data/upload/'
ALLOWED_EXTENSIONS = {'txt', 'tsv', 'csv'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/')
def index():
    return render_template('index.html')


@app.route('/documentation')
def documentation():
    return 'Documentation!  '


def get_config(form_results):
    parameters_cmd = "python3 pitviper_config.py "
    for parameter in form_results:
        value = form_results[parameter]
        if value == "":
            value = "none"
        parameters_cmd += " --{param} {value}".format(param=parameter, value=value)
    print(parameters_cmd)
    os.system(parameters_cmd)

def run_pitviper(token):
    configfile = 'config/{token}.yaml'.format(token=token)
    cmd = "python3 pitviper.py --run_snakemake True --configfile {conf}".format(conf=configfile)
    print(cmd)
    os.system(cmd)

@app.route('/result',methods = ['POST', 'GET'])
def result():
   if request.method == 'POST':
      f = request.files['library']
      print(">>>",f)
      f.save(secure_filename(f.filename))
      result = request.form
      # get_config(result)
      # run_pitviper(token=result['token'])
      return render_template("result.html",result = result)
