#!/home/paularthur/miniconda3/bin/python
import os
import subprocess
import webbrowser
from datetime import datetime
from threading import Timer

import yaml
from flask import Flask, render_template, request

app = Flask(__name__)


@app.route("/")
def index():
    """A Flask route that renders an HTML template file called "index.html"
    and returns it to the client."""
    return render_template("index.html")


@app.route("/documentation")
def documentation():
    """A Flask route that renders an HTML template file called "doc.html" and
    returns it to the client."""
    return render_template("doc.html")


def shutdown_server():
    """A function that shuts down the Flask development server."""
    func = request.environ.get("werkzeug.server.shutdown")
    if func is None:
        raise RuntimeError("Not running with the Werkzeug Server")
    func()


def run_pitviper(token):
    """A function that runs PitViper using a YAML configuration file
    and a token.
    Args:
    - token (str): A unique identifier for the run.

    Returns: None
    """
    configfile = f"config/{token}.yaml"
    with open(configfile, "r", encoding="utf-8") as stream:
        content = yaml.safe_load(stream)
    jobs = content.get("jobs", 1)
    print(f"Running PitViper with {jobs} jobs from the CLI...")
    cmd = f"python3 pitviper.py --configfile {configfile} --jobs {jobs} > logs/{token}/PitViper_snakemake.log 2>&1"
    print("Running command: ", cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as error:
        print("PitViper failed with error code: ", error.returncode)


@app.route("/result", methods=["POST", "GET"])
def result():
    """A Flask route that handles POST and GET requests and saves uploaded
    files to a resources directory.
    POST request is made when files are uploaded, GET request is made to
    retrieve the uploaded files.
    The function saves uploaded files to the resources directory and then
    writes a configuration YAML configuration file.
    Finally, it runs a function called run_pitviper with the token as an
    argument and shuts down the server.
    """
    result_dict = request.form.to_dict()

    def is_token_used(token):
        return os.path.exists(f"config/{token}.yaml")

    # Check if token isn't already used. If yes, add a timestamp to the token.
    token = result_dict["token"]
    token_used = is_token_used(token)
    while token_used:
        token_suffix = datetime.now().strftime("%d_%m_%Y_%H%M%S")
        new_token = f"{token}_{token_suffix}"
        if not os.path.exists(f"config/{new_token}.yaml"):
            token = new_token
            result_dict["token"] = token
            token_used = False
    if not os.path.exists(f"resources/{token}/"):
        os.makedirs(f"resources/{token}/")

    # Create logs directory if it does not exist
    if not os.path.exists(f"logs/{token}/"):
        os.makedirs(f"logs/{token}/")

    # Get uploaded files from POST request
    if request.method == "POST":
        library_file = request.files.get("library_file")
        controls_file = request.files.get("controls_file")
        bed_anno_file = request.files.get("bed_anno_file")
        tsv_file = request.files.get("tsv_file")
        count_table_file = request.files.get("count_table_file")

    # Create resources directory if it does not exist
    os.makedirs(f"resources/{result_dict['token']}/", exist_ok=True)

    # Save uploaded files in resources directory
    def save_file(filename, file):
        if file:
            filepath = f"resources/{result_dict['token']}/{filename}"
            file.save(filepath)
            return filepath
        else:
            return ""

    # Save uploaded files in resources directory
    result_dict["library_file"] = save_file(library_file.filename, library_file)
    result_dict["controls_file"] = save_file(controls_file.filename, controls_file)
    result_dict["bed_annotation_file"] = save_file("annotation.bed", bed_anno_file)
    result_dict["tsv_file"] = save_file(tsv_file.filename, tsv_file)

    # Save count table file in resources directory
    if count_table_file.filename != "":
        count_table_filename = (
            f"resources/{result_dict['token']}/{count_table_file.filename}"
        )
        count_table_file.save(count_table_filename)
        result_dict["count_table_file"] = count_table_filename
        result_dict[
            "normalized_count_table"
        ] = f"resources/{result_dict['token']}/screen.count_normalized.txt"
    else:
        result_dict[
            "count_table_file"
        ] = f"resources/{result_dict['token']}/screen.count.txt"
        result_dict[
            "normalized_count_table"
        ] = f"resources/{result_dict['token']}/screen.count_normalized.txt"

    # Save BAGEL files in resources directory if BAGEL is activated
    if result_dict["bagel_activate"] == "True":
        for filename in ["essentials", "nonessentials"]:
            file = request.files.get(filename)
            result_dict[filename] = save_file(file.filename, file)
    else:
        result_dict["essentials"] = ""
        result_dict["nonessentials"] = ""

    # Set mageck_count_activate based on bowtie_activate
    result_dict["mageck_count_activate"] = str(result_dict["bowtie_activate"] != "True")

    # Write YAML configuration file
    yaml_file_name = f"config/{result_dict['token']}.yaml"
    with open(yaml_file_name, "w", encoding="utf-8") as file:
        yaml.dump(result_dict, file)

    # Run PitViper and shutdown server, if PitViper fails, show error page
    try:
        run_pitviper(token=result_dict["token"])
        shutdown_server()
        if os.path.exists(f"results/{token}/Report.ipynb"):
            subprocess.check_call(f"jupyter notebook results/{token}", shell=True)
            return "You can close this window now."
        else:
            with open(
                f"logs/{token}/PitViper_snakemake.log", "r", encoding="utf-8"
            ) as f:
                log_content = f.read()
            return render_template(
                "error.html",
                error_message=f"PitViper failed to create results/{token}/.\
                Check logs/{token}/PitViper_snakemake.log for more information.",
                traceback=log_content,
            )
    except Exception as error:
        error_message = str(error)
        return render_template(
            "error.html",
            error_message=error_message,
            traceback="Error occured before PitViper was run.",
        )


def open_browser():
    """A function that opens the system's default web browser
    and navigates to http://127.0.0.1:5000/.
    """
    webbrowser.open_new("http://127.0.0.1:5000/")


# Start Flask server
Timer(1, open_browser).start()
