#!/home/paularthur/miniconda3/bin/python
import os
import webbrowser
from threading import Timer

import yaml
from flask import Flask, render_template, request

app = Flask(__name__)


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/documentation")
def documentation():
    return render_template("doc.html")


def shutdown_server():
    func = request.environ.get("werkzeug.server.shutdown")
    if func is None:
        raise RuntimeError("Not running with the Werkzeug Server")
    func()


def run_pitviper(token):
    configfile = f"config/{token}.yaml"
    with open(configfile, "r") as stream:
        content = yaml.safe_load(stream)
    print(content)
    cmd = f"python3 pitviper.py --configfile {configfile} --jobs {content['jobs']}"
    print(cmd)
    os.system(cmd)


@app.route("/result", methods=["POST", "GET"])
def result():
    print("PitViper is running from GUI.")
    result = request.form
    result_dict = dict(result)
    if request.method == "POST":
        try:
            library_file = request.files.get("library_file")
        except KeyError:
            library_file = False
        try:
            controls_file = request.files.get("controls_file")
        except KeyError:
            controls_file = False
        try:
            bed_anno_file = request.files.get("bed_anno_file")
        except KeyError:
            bed_anno_file = False

        tsv_file = request.files.get("tsv_file")
        count_table_file = request.files.get("count_table_file")

    # Check if resources/ directory exist, if not, create it.
    if not os.path.exists("resources/" + result_dict["token"] + "/"):
        os.makedirs("resources/" + result_dict["token"] + "/")

    # If a library file was uploaded, save it in resources/ directory.
    if library_file:
        library_filename = (
            "resources/" + result_dict["token"] + "/" + library_file.filename
        )
        library_file.save(library_filename)
    else:
        library_filename = ""
    result_dict["library_file"] = library_filename

    # If a control guides file was uploaded, save it in resources/ directory.
    if controls_file:
        controls_filename = (
            "resources/" + result_dict["token"] + "/" + controls_file.filename
        )
        controls_file.save(controls_filename)
    else:
        controls_filename = ""
    result_dict["controls_file"] = controls_filename

    # If an annotation file was uploaded, save it in resources/ directory.
    if bed_anno_file:
        bed_anno_filename = "resources/" + result_dict["token"] + "/annotation.bed"
        bed_anno_file.save(bed_anno_filename)
    else:
        bed_anno_filename = ""
    result_dict["bed_annotation_file"] = bed_anno_filename

    # Save design file in resources/ directory
    tsv_filename = "resources/" + result_dict["token"] + "/" + tsv_file.filename
    result_dict["tsv_file"] = tsv_filename
    tsv_file.save(tsv_filename)

    # Check if a count table was uploaded, if yes, save it in resources/ directory.
    if count_table_file.filename != "":
        count_table_filename = (
            "resources/" + result_dict["token"] + "/" + count_table_file.filename
        )
        count_table_file.save(count_table_filename)
        result_dict["count_table_file"] = count_table_filename
        result_dict["normalized_count_table"] = (
            "resources/" + result_dict["token"] + "/screen.count_normalized.txt"
        )
    else:
        result_dict["count_table_file"] = (
            "resources/" + result_dict["token"] + "/screen.count.txt"
        )
        result_dict["normalized_count_table"] = (
            "resources/" + result_dict["token"] + "/screen.count_normalized.txt"
        )

    # Check if BAGEL is activated, if not, create empty fields, otherwise save them in resources/ directory.
    if result_dict["bagel_activate"] == "False":
        result_dict["nonessentials"] = ""
        result_dict["essentials"] = ""
    else:
        ess = request.files["essentials"]
        ess_filename = "resources/" + result_dict["token"] + "/" + ess.filename
        ess.save(ess_filename)

        non_ess = request.files["nonessentials"]
        noness_filename = "resources/" + result_dict["token"] + "/" + non_ess.filename
        non_ess.save(noness_filename)

        result_dict["nonessentials"] = noness_filename
        result_dict["essentials"] = ess_filename

    # Init MAGeCK activation based on Bowtie2 value
    if result_dict["bowtie_activate"] == "True":
        result_dict["mageck_count_activate"] = "False"
    else:
        result_dict["mageck_count_activate"] = "True"

    # Write configuration YAML configuration file.
    yaml_file_name = "config/" + result["token"] + ".yaml"
    with open(yaml_file_name, "w") as file:
        documents = yaml.dump(result_dict, file)
    run_pitviper(token=result_dict["token"])
    # if success:
    #    Report file created
    shutdown_server()
    #  else:
    #     display failed page
    return "You can close this page and start using the Jupyter Notebook report."


def open_browser():
    webbrowser.open_new("http://127.0.0.1:5000/")


Timer(1, open_browser).start()
