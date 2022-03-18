#!/home/paularthur/miniconda3/bin/python
import os
import webbrowser
from threading import Timer

import yaml
from flask import Flask, redirect, render_template, request, url_for
from werkzeug.utils import secure_filename

app = Flask(__name__)


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/documentation")
def documentation():
    return "Documentation!"


def shutdown_server():
    func = request.environ.get("werkzeug.server.shutdown")
    if func is None:
        raise RuntimeError("Not running with the Werkzeug Server")
    func()


def run_pitviper(token):
    configfile = "config/{token}.yaml".format(token=token)
    with open(configfile, "r") as stream:
        content = yaml.safe_load(stream)
    print(content)
    cmd = "python3 pitviper.py --configfile {conf} --jobs {n}".format(
        conf=configfile, n=content["jobs"]
    )
    print(cmd)
    os.system(cmd)


@app.route("/result", methods=["POST", "GET"])
def result():
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

        tsv_file = request.files.get("tsv_file")
        count_table_file = request.files.get("count_table_file")

    if not os.path.exists("resources/" + result_dict["token"] + "/"):
        os.makedirs("resources/" + result_dict["token"] + "/")

    if library_file:
        library_filename = (
            "resources/" + result_dict["token"] + "/" + library_file.filename
        )
        library_file.save(library_filename)
    else:
        library_filename = ""
    result_dict["library_file"] = library_filename

    if controls_file:
        controls_filename = (
            "resources/" + result_dict["token"] + "/" + controls_file.filename
        )
        controls_file.save(controls_filename)
    else:
        controls_filename = ""
    result_dict["controls_file"] = controls_filename

    tsv_filename = "resources/" + result_dict["token"] + "/" + tsv_file.filename
    tsv_file.save(tsv_filename)

    result_dict["tsv_file"] = tsv_filename

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

    yaml_file_name = "config/" + result["token"] + ".yaml"
    with open(yaml_file_name, "w") as file:
        documents = yaml.dump(result_dict, file)
    run_pitviper(token=result_dict["token"])
    shutdown_server()
    return "You can close this page and start using the Jupyter Notebook report."


def open_browser():
    webbrowser.open_new("http://127.0.0.1:5000/")


Timer(1, open_browser).start()
