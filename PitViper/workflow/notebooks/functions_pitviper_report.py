import json
import os
import random as rd
import re
import warnings
import uuid
from functools import partial, reduce
from os import listdir, path
from pathlib import Path

import altair as alt
import IPython
import ipywidgets as widgets
import natsort as ns
import numpy as np
import pandas as pd
import plotly.figure_factory as ff
import plotly.graph_objects as go
import requests
import rpy2
import rpy2.ipython.html
import rpy2.robjects as ro
import yaml
from IPython.core.display import HTML, display
from rpy2.rinterface import RRuntimeWarning
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from rpy2.rinterface_lib.embedded import RRuntimeError
from scipy import stats
from scipy.stats import zscore
from sklearn import decomposition

from IPython.display import Markdown as md

import upsetplot
from upsetplot import from_memberships

import matplotlib.pyplot as plt


buf = []


def f(x):
    buf.append(x)


rpy2.rinterface_lib.callbacks.consolewrite_print = f
rpy2.rinterface_lib.callbacks.consolewrite_warnerror = f
rpy2.rinterface_lib.callbacks.showmessage = f

warnings.filterwarnings("ignore")

warnings.filterwarnings("ignore", category=RRuntimeWarning)

depmap = importr("depmap")
experimentHub = importr("ExperimentHub")
utils = importr("utils")

# Remove Altair max rows
alt.data_transformers.disable_max_rows()

pd.options.mode.chained_assignment = None

# Define layout for widgets
layout = widgets.Layout(width="auto", height="40px")  # set width and height
style = {"description_width": "initial"}


def display_warning(message):
    """Display a warning message in the notebook."""
    warning_html = """
    <div style='background-color: lightyellow; 
                border-left: 5px solid yellow; 
                padding: 10px; 
                margin-bottom: 10px;'>
        <p style='color: #f5aa42; 
                  font-weight: bold;'>{}</p>
    </div>
    """.format(
        message
    )
    display(HTML(warning_html))


def display_error(message):
    """Display an error message in the notebook."""
    error_html = """
    <div style='background-color: #ffdddd; 
                border-left: 5px solid red; 
                padding: 10px; 
                margin-bottom: 10px;'>
        <p style='color: red; 
                  font-weight: bold;'>{}</p>
    </div>
    """.format(
        message
    )
    display(HTML(error_html))


def display_info(message, bold=True):
    """Display an information message in the notebook."""
    font_weight = "bold" if bold else "normal"
    info_html = """
    <div style='background-color: #d0e7ff; 
                border-left: 5px solid blue; 
                padding: 10px; 
                margin: -0.4em;'>
        <p style='color: blue; 
                  font-weight: {};'>{}</p>
    </div>
    """.format(
        font_weight, message
    )
    display(HTML(info_html))


def natural_sort(l: list):
    """Do a natural sorting on the input list l.

    Args:
        l (list): List of strings.

    Returns:
        _type_: Natural sorted list of strings.
    """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
    return sorted(l, key=alphanum_key)


def working_directory_update(output: str):
    """Change working directory to ./PitViper

    Args:
        output (str): Name of the file to write.
    """
    switch = True
    while os.path.basename(os.getcwd()) != "PitViper":
        if switch:
            switch = False
            os.system("cd ../../")
        else:
            os.system("cd ../")


def open_yaml(yml: str):
    """Open and read content of YAML file yml.

    Args:
        yml (str): Path to yaml file to open and read.

    Returns:
        content (dict): Dictionnary containing yml content.
    """
    with open(yml, "r") as stream:
        try:
            content = yaml.safe_load(stream)
            return content
        except yaml.YAMLError as exc:
            print(exc)


def download(tools_available: dict, tool: str, treatment: str, control: str):
    """Display a button to download a file.

    Args:
        tools_available (dict): Dictionnary of PitViper results.
        tool (str): Tool name.
        treatment (str): Treatment condition.
        control (str): Control condition.
    """
    condition = f"{treatment}_vs_{control}"
    files_list = list(tools_available[tool][condition].keys())
    for file in files_list:
        content = tools_available[tool][condition][file].to_string()
        download_file(
            content=content, filename=tool + "_" + file, label=f"Download {file}"
        )


# def display_config(token: str):
#     config_name = "./config/%s.yaml" % token
#     config = open_yaml(config_name)
#     download_file(content=config, filename=config_name, label="Download config file!")


def download_config(token: str):
    config_name = f"./config/{token}.yaml"
    config_dict = open_yaml(config_name)

    if config_dict is not None:
        # Convert dictionary back to YAML formatted string
        config_yaml_str = yaml.dump(config_dict, default_flow_style=False)

        # Call download_file with the YAML string
        download_file(
            content=config_yaml_str, filename=config_name, label="Download config file!"
        )
    else:
        display_error("Error: Config file could not be loaded.")


def download_raw_counts(token):
    config_name = f"./config/{token}.yaml"
    config = open_yaml(config_name)

    # Get the name of the file with raw counts
    raw_cts_name = config["count_table_file"]

    # Load the file into a pandas DataFrame
    raw_cts_df = pd.read_table(raw_cts_name)

    # Convert DataFrame to TSV format string
    raw_cts_tsv = raw_cts_df.to_csv(sep="\t", index=False)

    # Call download_file with the TSV string
    download_file(
        content=raw_cts_tsv,
        filename=f"{raw_cts_name}.tsv",
        label="Download raw counts matrix!",
    )


def download_normalized_counts(token):
    config_name = f"./config/{token}.yaml"
    config = open_yaml(config_name)

    # Get the name of the normalized count table file
    cts_name = config["normalized_count_table"]

    # Load the file into a pandas DataFrame
    cts_df = pd.read_table(cts_name)

    # Convert DataFrame to TSV format string
    cts_tsv = cts_df.to_csv(sep="\t", index=False)

    # Call download_file with the TSV string
    download_file(
        content=cts_tsv,
        filename=f"{cts_name}.tsv",
        label="Download normalized counts matrix!",
    )


def download_design(token):
    config_name = "./config/%s.yaml" % token
    config = open_yaml(config_name)
    design_name = config["tsv_file"]
    design_content = pd.read_table(design_name).to_string()
    download_file(
        content=design_content,
        filename=design_name,
        label=f"Download design file!",
    )


def import_results(token: str, bagel_version: int = 2):
    """Load PitViper results inside token sub-directory.

    Args:
        token (str): Results token.

    Returns:
        Tuple(str, dict): Tuple containing directory scanned and dictionnary of results by tool, condition and file.
    """
    print("Token: %s\n" % token)
    config = "./config/%s.yaml" % token
    print("Config file used: %s" % config)
    content = open_yaml(config)

    tools = ["DESeq2"]

    if content["mageck_mle_activate"] == "True":
        tools.append("MAGeCK_MLE")
        print("\nPlease cite the following article if you use MAGeCK MLE:")
        print(
            "Li, W., Köster, J., Xu, H. et al. Quality control, modeling, and visualization of CRISPR screens with MAGeCK-VISPR. Genome Biol 16, 281 (2015). https://doi.org/10.1186/s13059-015-0843-6"
        )
    if content["mageck_rra_activate"] == "True":
        tools.append("MAGeCK_RRA")
        print("\nPlease cite the following article if you use MAGeCK RRA:")
        print(
            "Li, W., Xu, H., Xiao, T. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biol 15, 554 (2014). https://doi.org/10.1186/s13059-014-0554-4"
        )
    if content["bagel_activate"] == "True":
        if bagel_version != 1 and bagel_version != 2:
            raise ValueError("Invalid BAGEL version. Please use 1 or 2.")
        tools.append("BAGEL")
        print("\nPlease cite the following article if you use BAGEL:")
        print(
            "Li, W., Xu, H., Xiao, T. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biol 15, 554 (2014). https://doi.org/10.1186/s13059-014-0554-4"
        )
    if content["crisphiermix_activate"] == "True":
        tools.append("CRISPhieRmix")
        print("\nPlease cite the following article if you use CRISPhieRmix:")
        print(
            "Daley, T., Lin, Z., Lin, X. et al. CRISPhieRmix: a hierarchical mixture model for CRISPR pooled screens. Genome Biol 19, 159 (2018). https://doi.org/10.1186/s13059-018-1538-6"
        )
    if content["directional_scoring_method_activate"] == "True":
        tools.append("directional_scoring_method")
    if content["ssrea_activate"] == "True":
        tools.append("SSREA")

    results_directory = "results/%s/" % token
    print("\nResults directory: %s \n" % results_directory)

    tools_available = {}
    print("Tools available:")
    for tool in tools:
        if tool in os.listdir(results_directory):
            print("\t-%s" % tool)
            tools_available[tool] = {}

    for tool in tools_available:
        print("- Process %s results..." % tool)
        for comparison in os.listdir(results_directory + tool):
            tools_available[tool][comparison] = {}
            for file in os.listdir(os.path.join(results_directory, tool, comparison)):
                if file.endswith(".txt") or file.endswith(".pr") or file.endswith(".bf"):
                    if tool in ["CRISPhieRmix"]:
                        sep = ","
                    else:
                        sep = "\t"
                    tools_available[tool][comparison][file] = pd.read_csv(
                        os.path.join(results_directory, tool, comparison, file), sep=sep
                    )
    add_columns(tools_available)
    return (results_directory, tools_available)


def add_columns(tools_available: dict):
    """Add columns to result dataframes in tools_available.

    Args:
        tools_available (dict): Dictionnary of PitViper results.
    """
    if "MAGeCK_MLE" in tools_available:
        for condition in tools_available["MAGeCK_MLE"]:
            treatment = condition.split("_vs_")[0]
            for file_suffixe in ["%s.gene_summary.txt", "%s.sgrna_summary.txt"]:
                table = tools_available["MAGeCK_MLE"][condition][
                    file_suffixe % condition
                ]
                if (not "log10(invFDR)" in list(table.columns)) and (
                    "%s|fdr" % treatment in list(table.columns)
                ):
                    array = table["%s|fdr" % treatment].values
                    min_fdr = np.min(array[np.nonzero(array)])
                    table["%s|fdr_nozero" % treatment] = table[
                        "%s|fdr" % treatment
                    ].replace(0, min_fdr)
                    table["log10(invFDR)"] = -np.log10(
                        table["%s|fdr_nozero" % treatment]
                    )
    if "MAGeCK_RRA" in tools_available:
        for condition in tools_available["MAGeCK_RRA"]:
            treatment = condition.split("_vs_")[0]
            for file_suffixe in ["%s.gene_summary.txt", "%s.sgrna_summary.txt"]:
                table = tools_available["MAGeCK_RRA"][condition][
                    file_suffixe % condition
                ]
                for direction in ["neg", "pos"]:
                    if not "%s|log10(invFDR)" % direction in list(table.columns) and (
                        "%s|fdr" % direction in list(table.columns)
                    ):
                        array = table["%s|fdr" % direction].values
                        min_fdr = np.min(array[np.nonzero(array)])
                        table["%s|fdr_nozero" % direction] = table[
                            "%s|fdr" % direction
                        ].replace(0, min_fdr)
                        table["%s|log10(invFDR)" % direction] = -np.log10(
                            table["%s|fdr_nozero" % direction]
                        )
    if "SSREA" in tools_available:
        for condition in tools_available["SSREA"]:
            treatment = condition.split("_vs_")[0]
            for file_suffixe in ["%s_all-elements_SSREA.txt"]:
                table = tools_available["SSREA"][condition][file_suffixe % condition]
                if (not "log10(invPadj)" in list(table.columns)) and (
                    "padj" in list(table.columns)
                ):
                    array = table["padj"].values
                    min_fdr = np.min(array[np.nonzero(array)])
                    table["padj_nozero"] = table["padj"].replace(0, min_fdr)
                    table["log10(invPadj)"] = -np.log10(table["padj_nozero"])
    if "CRISPhieRmix" in tools_available:
        for condition in tools_available["CRISPhieRmix"]:
            treatment = condition.split("_vs_")[0]
            for file_suffixe in ["%s.txt"]:
                table = tools_available["CRISPhieRmix"][condition][
                    file_suffixe % condition
                ]
                if not "log10(invFDR)" in list(table.columns):
                    array = table["locfdr"].values
                    min_fdr = np.min(array[np.nonzero(array)])
                    table["padj_nozero"] = table["locfdr"].replace(0, min_fdr)
                    table["log10(invPadj)"] = -np.log10(table["padj_nozero"])


def show_mapping_qc(token: str):
    """Display mapping quality control table, if available.

    Args:
        token (str): Results token.

    Returns:
        pandas.DataFrame: DataFrame of mapping quality control.
    """
    path_qc = "./resources/%s/screen.countsummary.txt" % token
    if not path.exists(path_qc):
        display_warning("No mapping QC file to show.")
        # print("No mapping QC file to show.")
        return 0
    table = pd.read_csv(path_qc, sep="\t")
    table = table[["Label", "Reads", "Mapped", "Percentage", "Zerocounts", "GiniIndex"]]
    table["Label"] = pd.Categorical(
        table["Label"], ordered=True, categories=ns.natsorted(table["Label"].unique())
    )
    table = table.sort_values("Label")

    def color_low_mapping_red(val):
        color = "red" if float(val) < 0.6 else "green"
        return "color: %s" % color

    def color_high_gini_red(val):
        color = "red" if val > 0.35 else "green"
        return "color: %s" % color

    s = table.style.applymap(color_low_mapping_red, subset=["Percentage"]).applymap(
        color_high_gini_red, subset=["GiniIndex"]
    )
    display(s)


# def download_file(content, label, filename):
#     # Add download button
#     outname = os.path.basename(filename)
#     id_file = rd.random()
#     display(
#         HTML(
#             '<textarea id="textbox_{id_file}" style="display: none;">{content}</textarea> <button id="create_{id_file}">{label}</button> <a download="{filename}" id="downloadlink_{id_file}" style="display: none">Download</a>'.format(
#                 **locals()
#             )
#         )
#     )
#     display(
#         HTML(
#             '<script type="text/javascript">!function(){{var e=null,t=document.getElementById("create_{id_file}"),n=document.getElementById("textbox_{id_file}");t.addEventListener("click",function(){{var t,l,c=document.getElementById("downloadlink_{id_file}");c.href=(t=n.value,l=new Blob([t],{{type:"text/plain"}}),null!==e&&window.URL.revokeObjectURL(e),e=window.URL.createObjectURL(l)),c.click()}},!1)}}();</script>'.format(
#                 **locals()
#             )
#         )
#     )


def download_file(content, label, filename):
    """
    Generates a download button in a Jupyter Notebook for the given content.

    Parameters:
    content (str): Content to be downloaded.
    label (str): Label for the download button.
    filename (str): Name of the file to be downloaded.

    Returns:
    None
    """
    try:
        # Generate a unique ID for HTML elements
        unique_id = str(uuid.uuid4())

        # Get the base name of the file to ensure it's a valid filename
        safe_filename = os.path.basename(filename)

        # Determine the MIME type (default to 'text/plain')
        mime_type = "text/plain"

        # HTML and JavaScript code for the download button
        html_script = f"""
        <textarea id="textbox_{unique_id}" style="display: none;">{content}</textarea>
        <button id="create_{unique_id}">{label}</button>
        <a download="{safe_filename}" id="downloadlink_{unique_id}" style="display: none">Download</a>
        <script type="text/javascript">
        !function() {{
            var file_blob = null;
            var button = document.getElementById("create_{unique_id}");
            var textbox = document.getElementById("textbox_{unique_id}");
            button.addEventListener("click", function() {{
                var content = textbox.value;
                var blob = new Blob([content], {{type: "{mime_type}"}});
                if (file_blob !== null) {{
                    window.URL.revokeObjectURL(file_blob);
                }}
                file_blob = window.URL.createObjectURL(blob);
                var downloadLink = document.getElementById("downloadlink_{unique_id}");
                downloadLink.href = file_blob;
                downloadLink.click();
            }}, false);
        }}();
        </script>
        """

        # Display the HTML and JavaScript
        display(HTML(html_script))
    except Exception as e:
        display_error(f"An error occurred: {e}")
        # print(f"An error occurred: {e}")


def show_read_count_distribution(token: str, width=800, height=400):
    """Display an altair chart of read counts distribution.

    Args:
        token (str): Results token.
        width (int, optional): Figure width. Defaults to 800.
        height (int, optional): Figure height. Defaults to 400.

    Returns:
        altair.Chart: Chart of normalized read counts distribution.
    """
    config = "./config/%s.yaml" % token
    content = open_yaml(config)

    # Load the design file
    design = pd.read_csv(content["tsv_file"], sep="\t")

    # Create a dictionary to map replicate names to condition names
    d = dict(zip(design.replicate, design.condition))

    # Load the normalized count table
    norm_file = content["normalized_count_table"]
    if not path.exists(norm_file):
        display_info("No count file to show.")
        # print("No count file to show.")
        return 0
    table = pd.read_csv(norm_file, sep="\t")

    table = table[[col for col in d]]

    table += 1
    table = table.apply(np.log2)

    chart = (
        (
            alt.Chart(table)
            .transform_fold(list(table.columns), as_=["Replicate", "counts"])
            .transform_density(
                density="counts",
                bandwidth=0.3,
                groupby=["Replicate"],
                extent=[0, 20],
                counts=True,
                steps=200,
            )
            .mark_line()
            .encode(
                alt.X("value:Q", axis=alt.Axis(title="log2(read count)")),
                alt.Y("density:Q", axis=alt.Axis(title="Density")),
                alt.Color("Replicate:N"),
                tooltip=["Replicate:N", "value:Q", "density:Q"],
            )
            .properties(width=width, height=height)
        )
        .configure_axis(grid=False)
        .configure_view(stroke=None)
    )

    return chart


def pca_counts(token: str):
    """Compute and display PCA of normalized read counts.

    Args:
        token (str): Results token.

    Returns:
        altair.Chart: Altair chart
    """
    config = "./config/%s.yaml" % token
    content = open_yaml(config)

    # Load the design file
    design = pd.read_csv(content["tsv_file"], sep="\t")

    # Create a dictionary to map replicate names to condition names
    d = dict(zip(design.replicate, design.condition))

    # Load the normalized count table
    cts_file = content["normalized_count_table"]
    cts = pd.read_csv(cts_file, sep="\t")

    # Filter columns to keep only replicates found in the design file
    cts = cts[[col for col in d]]

    # Create a numpy array of condition names
    y = np.array([d[k] for k in cts.columns if k in d])
    # Create a numpy array of replicate names
    y_bis = np.array(cts.columns)

    # Transpose the count table
    X = cts.to_numpy().T

    # Perform PCA
    pca = decomposition.PCA(n_components=2)
    pca.fit(X)
    X = pca.transform(X)

    # Create a pandas DataFrame from the PCA results
    a = pd.DataFrame(X, columns=["PC1", "PC2"])
    b = pd.DataFrame(y, columns=["condition"])
    c = pd.DataFrame(y_bis, columns=["replicate"])

    source = pd.concat([a, b, c], axis=1)

    PC1_explained_variance_ratio = round(pca.explained_variance_ratio_[0] * 100, 2)
    PC2_explained_variance_ratio = round(pca.explained_variance_ratio_[1] * 100, 2)

    pca_2d = (
        (
            alt.Chart(source)
            .mark_circle(size=60)
            .encode(
                x=alt.X(
                    "PC1:Q",
                    axis=alt.Axis(
                        title="PC1 ({p}%)".format(p=PC1_explained_variance_ratio)
                    ),
                ),
                y=alt.X(
                    "PC2:Q",
                    axis=alt.Axis(
                        title="PC2 ({p}%)".format(p=PC2_explained_variance_ratio)
                    ),
                ),
                color="condition:N",
                tooltip=["PC1", "PC2", "condition", "replicate"],
            )
            .interactive()
        )
        .configure_axis(grid=False)
        .configure_view(stroke=None)
    )

    return pca_2d


def get_enrichr_results(genes: list, description: str, gene_set_library: list):
    """Get and return EnrichR results for a given list of genes.

    Args:
        genes (list): List of genes to use for EnrichR analysis.
        description (str): Description of the genes list.
        gene_set_library (list): List of genes set to use.

    Returns:
        dict: EnrichR results, or None in case of errors
    """
    enrichr_url_addlist = "http://maayanlab.cloud/Enrichr/addList"
    genes_str = "\n".join(genes)
    payload = {"list": (None, genes_str), "description": (None, description)}

    try:
        # Send POST request to EnrichR API
        response = requests.post(enrichr_url_addlist, files=payload, timeout=10)
        response.raise_for_status()  # Will raise an exception for HTTP errors
    except requests.exceptions.RequestException as e:
        # Print error message and return None
        display_error(f"Error analyzing gene list: {e}")
        # print("Error analyzing gene list: %s", e)
        return None

    try:
        # Parse JSON response
        data = response.json()
    except json.JSONDecodeError:
        display_error("Invalid JSON response")
        # print("Invalid JSON response")
        return None

    if "userListId" not in data:
        # Print error message and return None
        # print("Missing 'userListId' in response")
        display_error("Missing 'userListId' in response")
        return None

    user_list_id = data["userListId"]

    enrichr_url_addlist = "http://maayanlab.cloud/Enrichr/enrich"
    query_string = "?userListId=%s&backgroundType=%s"
    try:
        response = requests.get(
            enrichr_url_addlist + query_string % (user_list_id, gene_set_library),
            timeout=10,
        )
        response.raise_for_status()  # Will raise an exception for HTTP errors
    except requests.exceptions.RequestException as e:
        # print("Error fetching enrichment results: %s", e)
        display_error(f"Error fetching enrichment results: {e}")
        return None

    try:
        data = response.json()
    except json.JSONDecodeError:
        # print("Invalid JSON response")
        display_error("Invalid JSON response")
        return None

    return data


def create_enrichr_table(enrichrResults: dict):
    """Convert enrichr results dict to a pandas DataFrame.

    Args:
        enrichrResults (dict): EnrichR results

    Returns:
        pandas.DataFrame: EnrichR results as a pandas dataframe.
    """
    cols = [
        "Rank",
        "Term name",
        "P-value",
        "Z-score",
        "Combined score",
        "Overlapping genes",
        "Adjusted p-value",
        "Old p-value",
        "Old adjusted p-value",
    ]
    for base in enrichrResults:
        rows = []
        for pathway in enrichrResults[base]:
            row = pd.DataFrame([pathway])
            rows.append(row)

    table = pd.concat(rows)
    table.columns = cols
    return table


def enrichmentBarPlot(source, n, description, col_1, col_2, base):
    if n == "max":
        n = len(source.index)
    source = source.sort_values(by=["Combined score"], ascending=False).head(n)

    domain = [source["Adjusted p-value"].min(), source["Adjusted p-value"].max()]
    range_ = [col_1, col_2]

    bars = (
        alt.Chart(source)
        .mark_bar()
        .encode(
            x="Combined score",
            y=alt.Y("Term name", sort="-x"),
            tooltip=[
                "Rank",
                "Term name",
                "P-value",
                "Z-score",
                "Combined score",
                "Overlapping genes",
                "Adjusted p-value",
                "Old p-value",
                "Old adjusted p-value",
            ],
            color=alt.Color(
                "Adjusted p-value",
                scale=alt.Scale(domain=domain, range=range_),
                legend=alt.Legend(title="Adjusted p-value:"),
            ),
        )
        .properties(
            title=f"{description} ({base})",
        )
    )

    chart = (bars).properties(height=15 * n, width=500)
    return chart


def enrichmentCirclePlot(source, n, description, col_1, col_2, base):
    if n == "max":
        n = int(len(source.index))

    source = source.sort_values(by=["Combined score"], ascending=False).head(n)

    source["Overlap size"] = source["Overlapping genes"].str.len()

    domain = [source["Adjusted p-value"].min(), source["Adjusted p-value"].max()]
    range_ = [col_1, col_2]

    chart = (
        alt.Chart(source)
        .mark_circle(size=60)
        .encode(
            alt.X(
                "Combined score",
                scale=alt.Scale(
                    domain=(
                        source.min()["Combined score"],
                        source.max()["Combined score"] * 1.05,
                    )
                ),
            ),
            y=alt.Y("Term name", sort="-x"),
            color=alt.Color(
                "Adjusted p-value",
                scale=alt.Scale(domain=domain, range=range_),
                legend=alt.Legend(title="Adjusted p-value:"),
            ),
            tooltip=[
                "Rank",
                "Term name",
                "P-value",
                "Z-score",
                "Combined score",
                "Overlapping genes",
                "Adjusted p-value",
                "Old p-value",
                "Old adjusted p-value",
                "Overlap size",
            ],
            size="Overlap size",
        )
        .properties(title=description + " (%s)" % base[:-1])
        .configure_axisY(titleAngle=0, titleY=-10, titleX=-60, labelLimit=250)
        .interactive()
    )

    chart = (chart).properties(height=20 * n, width=500)

    return chart


def SSREA_like_data(
    comparison="", control="", tool="", results_directory="", tools_available=""
):
    """Return SSREA results as pandas dataframe."""
    tables_list = []
    for _comparison in os.listdir(os.path.join(results_directory, tool)):
        if _comparison.split("_vs_")[-1] == control:
            keys_list = list(tools_available[tool][_comparison].keys())
            for key in keys_list:
                if key.endswith("_all-elements_SSREA.txt"):
                    break
            data = tools_available[tool][_comparison][key]
            trt = _comparison.split("_vs_")[0]
            data["condition"] = trt
            tables_list.append(data)
        if _comparison == comparison:
            data = tools_available[tool][_comparison][
                "%s_all-elements_SSREA.txt" % comparison
            ]
            tables_list.append(data)
    result = pd.concat(tables_list)
    return result


def directional_scoring_method_data(
    comparison="", control="", tool="", results_directory="", tools_available=""
):
    """Return Directional Scoring Method results as pandas dataframe."""
    tables_list = []
    for _comparison in os.listdir(os.path.join(results_directory, tool)):
        if _comparison.split("_vs_")[-1] == control:
            keys_list = list(tools_available[tool][_comparison].keys())
            for key in keys_list:
                if key.endswith("_all-elements_directional_scoring_method.txt"):
                    break
            data = tools_available[tool][_comparison][key]
            trt = _comparison.split("_vs_")[0]
            data["condition"] = trt
            tables_list.append(data)
        if _comparison == comparison:
            data = tools_available[tool][_comparison][
                "%s_all-elements_directional_scoring_method.txt" % comparison
            ]
            tables_list.append(data)
    result = pd.concat(tables_list)
    return result


def MAGeCK_RRA_data(
    comparison="",
    control="",
    tool="MAGeCK_RRA",
    results_directory="",
    tools_available="",
):
    """Return MAGeCK RRA results as pandas dataframe."""

    def check(comparison, _comparison, control, mode):
        if mode:
            if _comparison.split("_vs_")[-1] == control:
                return True
            else:
                return False
        else:
            if _comparison == comparison:
                return True
            else:
                return False

    tables_list = []
    if control != "" and comparison == "":
        mode = True
    elif comparison != "" and control == "":
        mode = False
    else:
        display_error("Error, please check your inputs.")
    for _comparison in os.listdir(os.path.join(results_directory, tool)):
        if check(comparison, _comparison, control, mode):
            keys_list = list(tools_available[tool][_comparison].keys())
            for key in keys_list:
                if key.endswith(".gene_summary.txt"):
                    break
            data = tools_available[tool][_comparison][key]
            trt = _comparison.split("_vs_")[0]
            data["condition"] = trt
            tables_list.append(data)
    result = pd.concat(tables_list)
    return result


def MAGeCK_MLE_data(
    comparison="", control="", tool="", results_directory="", tools_available=""
):
    """Return MAGeCK MLE results as pandas dataframe."""

    def check(comparison, _comparison, control, mode):
        if mode:
            if _comparison.split("_vs_")[-1] == control:
                return True
            else:
                return False
        else:
            if _comparison == comparison:
                return True
            else:
                return False

    tables_list = []
    if control != "" and comparison == "":
        mode = True
    elif comparison != "" and control == "":
        mode = False
    for _comparison in os.listdir(os.path.join(results_directory, tool)):
        if check(comparison, _comparison, control, mode):
            keys_list = list(tools_available[tool][_comparison].keys())
            for key in keys_list:
                if key.endswith(".gene_summary.txt"):
                    break
            data = tools_available[tool][_comparison][key]
            trt = _comparison.split("_vs_")[0]
            data["condition"] = trt
            if mode:
                data = data[data.columns.drop(list(data.filter(regex="%s" % control)))]
                data = data.rename(columns=lambda x: re.sub(".+\|", "", x))
                tables_list.append(data)
            else:
                tables_list.append(data)
    result = pd.concat(tables_list)
    return result


def BAGEL_data(
    comparison="", control="", tool="BAGEL", results_directory="", tools_available="", bagel_version=2,
):
    """Return BAGEL results as pandas dataframe."""
    def check(comparison, _comparison, control, mode):
        if mode:
            if _comparison.split("_vs_")[-1] == control:
                return True
            else:
                return False
        else:
            if _comparison == comparison:
                return True
            else:
                return False

    if bagel_version == 2:
        end_file = "_BAGEL_output.pr"
    elif bagel_version == 1:
        end_file = "_BAGEL1_output.bf"

    tables_list = []
    if control != "" and comparison == "":
        mode = True
    elif comparison != "" and control == "":
        mode = False
    else:
        display_error("Error, please check your inputs.")
    for _comparison in os.listdir(os.path.join(results_directory, tool)):
        if _comparison.split("_vs_")[-1] == control:
            keys_list = list(tools_available[tool][_comparison].keys())
            for key in keys_list:
                if key.endswith(end_file):
                    break
            data = tools_available[tool][_comparison][key]
            trt = _comparison.split("_vs_")[0]
            data["condition"] = trt
            tables_list.append(data)
        if _comparison == comparison:
            data = tools_available[tool][_comparison][
                _comparison + end_file
            ]
            tables_list.append(data)
    result = pd.concat(tables_list)
    return result


def CRISPhieRmix_data(
    comparison="", control="", tool="", results_directory="", tools_available=""
):
    """Return CRISPhieRmix results as pandas dataframe."""
    tables_list = []
    for _comparison in os.listdir(os.path.join(results_directory, tool)):
        if _comparison.split("_vs_")[-1] == control:
            keys_list = list(tools_available[tool][_comparison].keys())
            data = tools_available[tool][_comparison][keys_list[0]]
            trt = _comparison.split("_vs_")[0]
            data["condition"] = trt
            tables_list.append(data)
        if _comparison == comparison:
            data = tools_available[tool][_comparison]["%s.txt" % comparison]
            tables_list.append(data)
    result = pd.concat(tables_list)
    return result


def filter_by_threshold(
    df,
    score_column,
    threshold,
    orientation,
    significant_label,
    fdr_column=None,
    fdr_cutoff=None,
):
    """Filter a dataframe based on a threshold and a column."""
    orientation = {
        ">=": lambda x: x > threshold,
        "<=": lambda x: x <= threshold,
        "abs() >=": lambda x: abs(x) >= threshold,
        "abs() <=": lambda x: abs(x) <= threshold,
        "abs() >": lambda x: abs(x) > threshold,
        "abs() <": lambda x: abs(x) < threshold,
    }[orientation]

    if fdr_cutoff and fdr_column:
        mask = (df[fdr_column] < fdr_cutoff) & (orientation(df[score_column]))
    else:
        mask = orientation(df[score_column])

    df.loc[mask, "significant"] = significant_label
    return df


def get_reverse_orientation(orientation):
    """Return the reverse orientation of a given orientation."""
    if orientation == ">=":
        return "<="
    elif orientation == "<=":
        return ">="
    elif orientation == "abs() >=":
        return "abs() <"
    elif orientation == "abs() <=":
        return "abs() >"
    elif orientation == "abs() >":
        return "abs() <="
    elif orientation == "abs() <":
        return "abs() >="


def get_pretty_orientation(orientation):
    """Return the pretty version of a given orientation."""
    if orientation == ">=":
        return "≥"
    elif orientation == "<=":
        return "≤"
    elif orientation == "abs() >=":
        return "abs≥"
    elif orientation == "abs() <=":
        return "abs≤"
    elif orientation == "abs() >":
        return "abs() >"
    elif orientation == "abs() <":
        return "abs() <"


def tool_results(results_directory, tools_available, token, bagel_version):
    """Display selected method's results for all genes."""
    config = f"./config/{token}.yaml"
    content = open_yaml(config)
    cts_file = "results/%s/normalized.filtered.counts.txt" % token
    cts = pd.read_csv(cts_file, sep="\t")
    cts_columns = [col for col in cts.columns.tolist() if not col in ["sgRNA", "Gene"]]
    cts = pd.melt(cts, id_vars=["sgRNA", "Gene"])
    genes_list = list(set(cts.Gene.tolist()))

    tools = [tool for tool in tools_available.keys() if tool != "DESeq2"]
    tools_widget = widgets.SelectMultiple(
        options=set(tools), description="Tool:", value=(tools[0],)
    )

    def update_comparisons_widget(new):
        comparisons_widget.options = tools_available[tools_widget.value[0]].keys()

    tools_widget.observe(update_comparisons_widget, "value")

    comparisons_list = os.listdir(
        os.path.join(results_directory, tools_widget.value[0])
    )
    comparisons_widget = widgets.Dropdown(
        options=set(comparisons_list), description="Comparison:"
    )

    element = widgets.Combobox(
        placeholder="Choose one",
        options=genes_list,
        description="Element(s):",
        ensure_option=False,
    )

    fdr_widget = widgets.FloatSlider(
        min=0.0, max=1.0, step=0.01, value=0.05, description="FDR cut-off:"
    )

    # Add a widget to define the BAGEL minimum BF cut-off. Default is 0. No minimum and maximum values are defined. Using FloatText widget.
    bagel_bf_widget = widgets.FloatText(
        value=0,
        description="BAGEL BF cut-off:",
        # If BAGEL is not available, the widget is disabled
        disabled="BAGEL" not in tools,
        style=style,
    )

    # Add a widget to define the orientation of the BAGEL BF cut-off. Default is ">=".
    bagel_bf_orientation_widget = widgets.Dropdown(
        options=[">=", "<=", "abs() >=", "abs() <="],
        value=">=",
        description="",
        # If BAGEL is not available, the widget is disabled
        disabled="BAGEL" not in tools,
        # style=style,
        layout=widgets.Layout(width="75px"),
    )

    # Add a widget to define the SSREA minimum NES cut-off. Default is 0. No minimum and maximum values are defined. Using FloatText widget.
    ssrea_nes_widget = widgets.FloatText(
        value=0,
        description="SSREA NES cut-off:",
        # If SSREA is not available, the widget is disabled
        disabled="SSREA" not in tools,
        style=style,
    )

    # Add a widget to define the orientation of the SSREA NES cut-off. Default is ">=". Using Dropdown widget.
    ssrea_nes_orientation_widget = widgets.Dropdown(
        options=[">=", "<=", "abs() >=", "abs() <="],
        value="abs() >=",
        description="",
        # If SSREA is not available, the widget is disabled
        disabled="SSREA" not in tools,
        # style=style,
        layout=widgets.Layout(width="75px"),
    )

    # Add a widget to define the directional_scoring_method minimum score cut-off. Default is 0. No minimum and maximum values are defined. Using FloatText widget.
    directional_scoring_method_score_widget = widgets.FloatText(
        value=0,
        description="DSM score cut-off:",
        # If directional_scoring_method is not available, the widget is disabled
        disabled="directional_scoring_method" not in tools,
        style=style,
    )

    # Add a widget to define the orientation of the directional_scoring_method score cut-off. Default is ">=". Using Dropdown widget.
    directional_scoring_method_score_orientation_widget = widgets.Dropdown(
        options=[">=", "<=", "abs() >=", "abs() <=", "abs() >", "abs() <"],
        value="abs() >",
        description="",
        # If directional_scoring_method is not available, the widget is disabled
        disabled="directional_scoring_method" not in tools,
        # style=style,
        layout=widgets.Layout(width="75px"),
    )

    # Add a widget to define the CRISPhieRmix minimum mean logFC cut-off. Default is 0. No minimum and maximum values are defined. Using FloatText widget.
    crisphiermix_logfc_widget = widgets.FloatText(
        value=0,
        description="CRISPhieRmix logFC cut-off:",
        # If CRISPhieRmix is not available, the widget is disabled
        disabled="CRISPhieRmix" not in tools,
        style=style,
    )

    # Add a widget to define the orientation of the CRISPhieRmix mean logFC cut-off. Default is ">=". Using Dropdown widget.
    crisphiermix_logfc_orientation_widget = widgets.Dropdown(
        options=[">=", "<=", "abs() >=", "abs() <="],
        value="abs() >=",
        description="",
        # If CRISPhieRmix is not available, the widget is disabled
        disabled="CRISPhieRmix" not in tools,
        # style=style,
        layout=widgets.Layout(width="75px"),
    )

    # Add a widget to define the MAGeCK_MLE minimum beta cut-off. Default is 0. No minimum and maximum values are defined. Using FloatText widget.
    mageck_mle_beta_widget = widgets.FloatText(
        value=0,
        description="MAGeCK MLE beta cut-off:",
        # If MAGeCK_MLE is not available, the widget is disabled
        disabled="MAGeCK_MLE" not in tools,
        style=style,
    )

    # Add a widget to define the orientation of the MAGeCK_MLE beta cut-off. Default is ">=". Using Dropdown widget.
    mageck_mle_beta_orientation_widget = widgets.Dropdown(
        options=[">=", "<=", "abs() >=", "abs() <="],
        value="abs() >=",
        description="",
        # If MAGeCK_MLE is not available, the widget is disabled
        disabled="MAGeCK_MLE" not in tools,
        # style=style,
        layout=widgets.Layout(width="75px"),
    )

    # Add a widget to define the MAGeCK_RRA minimum LFC cut-off. Default is 0. No minimum and maximum values are defined. Using FloatText widget.
    mageck_rra_lfc_widget = widgets.FloatText(
        value=0,
        description="MAGeCK RRA LFC cut-off:",
        # If MAGeCK_RRA is not available, the widget is disabled
        disabled="MAGeCK_RRA" not in tools,
        style=style,
    )

    # Add a widget to define the orientation of the MAGeCK_RRA LFC cut-off. Default is ">=". Using Dropdown widget.
    mageck_rra_lfc_orientation_widget = widgets.Dropdown(
        options=[">=", "<=", "abs() >=", "abs() <="],
        value="abs() >=",
        description="",
        # If MAGeCK_RRA is not available, the widget is disabled
        disabled="MAGeCK_RRA" not in tools,
        # style=style,
        layout=widgets.Layout(width="75px"),
    )

    color_sig_widget = widgets.ColorPicker(
        concise=False,
        description="Significant color:",
        value="red",
        style=style,
    )
    color_non_widget = widgets.ColorPicker(
        concise=False,
        description="Non-significant color:",
        value="gray",
        style=style,
    )

    def _MAGeCK_MLE_snake_plot(
        comparison,
        fdr_cutoff,
        non_sig,
        sig,
        results_directory,
        tools_available,
        elements,
    ):
        # Define the tool
        tool = "MAGeCK_MLE"

        # Define the beta and beta orientation
        mageck_mle_beta = float(mageck_mle_beta_widget.value)
        mageck_mle_beta_orientation = mageck_mle_beta_orientation_widget.value

        # Define the significant and non-significant labels using FDR and beta thresholds
        significant_label = f"FDR < {fdr_cutoff} and beta {get_pretty_orientation(mageck_mle_beta_orientation)} {mageck_mle_beta}"

        # Define the non-significant label
        non_significant_label = f"FDR ≥ {fdr_cutoff} or beta {get_pretty_orientation(get_reverse_orientation(mageck_mle_beta_orientation))} {mageck_mle_beta}"

        # Define the highlight label
        highlight_label = "Hit(s) of Interest"

        # Define the treatment and control names
        treatment, control = comparison.split("_vs_")

        # Get the MAGeCK MLE results
        source = MAGeCK_MLE_data(
            comparison=comparison,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )

        # Compute the default rank of the beta scores
        source["default_rank"] = source[treatment + "|beta"].rank()

        # Set the significant label to non-significant by default
        source["significant"] = non_significant_label

        # Filter the data based on the beta and FDR thresholds
        source = filter_by_threshold(
            source,
            treatment + "|beta",
            mageck_mle_beta,
            mageck_mle_beta_orientation,
            significant_label,
            treatment + "|fdr",
            fdr_cutoff,
        )

        # Highlight the elements of interest
        source.loc[source.Gene.isin(elements), "significant"] = highlight_label

        # Define the domain and range for the color scale
        domain = [significant_label, non_significant_label, highlight_label]
        range_ = [sig, non_sig, "blue"]

        # Define the color scale
        source["significant"] = pd.Categorical(
            source["significant"],
            categories=[non_significant_label, significant_label, highlight_label],
            ordered=True,
        )

        # Order rows by label 'significant' to have the highlighted genes on top
        source = source.sort_values(by="significant", ascending=True)

        # Create the snake plot
        chart = (
            alt.Chart(source, title=f"MAGeCK MLE ({comparison})")
            .mark_circle(size=60)
            .encode(
                x=alt.X("default_rank:Q", axis=alt.Axis(title="Rank")),
                y=alt.Y(
                    treatment + "|beta:Q", axis=alt.Axis(title=f"{treatment} beta")
                ),
                tooltip=[
                    "Gene",
                    "sgRNA",
                    treatment + "|beta",
                    treatment + "|fdr",
                    treatment + "|p-value",
                    "significant",
                    "default_rank",
                ],
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
            )
            .properties(width=800, height=400)
            .interactive()
        )
        line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule().encode(y="y")
        text = (
            alt.Chart(source.query("significant == 'Hit(s) of Interest'"))
            .mark_text(dy=10, dx=20, color="blue")
            .encode(
                x=alt.X("default_rank:Q"),
                y=alt.Y(treatment + "|beta:Q"),
                text=alt.Text("Gene"),
            )
        )
        chart = (
            (chart + line + text).configure_axis(grid=False).configure_view(stroke=None)
        )
        display(chart)

    def _MAGeCK_RRA_snake_plot(
        comparison,
        fdr_cutoff,
        non_sig,
        sig,
        results_directory,
        tools_available,
        elements,
    ):
        # Define the tool
        tool = "MAGeCK_RRA"

        # Define the LFC and LFC orientation
        mageck_rra_lfc = float(mageck_rra_lfc_widget.value)
        mageck_rra_lfc_orientation = mageck_rra_lfc_orientation_widget.value

        # Define the significant and non-significant labels using FDR threshold and LFC
        significant_label = f"FDR < {fdr_cutoff} and LFC {get_pretty_orientation(mageck_rra_lfc_orientation)} {mageck_rra_lfc}"
        non_significant_label = f"FDR ≥ {fdr_cutoff} or LFC {get_pretty_orientation(get_reverse_orientation(mageck_rra_lfc_orientation))} {mageck_rra_lfc}"

        # Define the highlight label
        highlight_label = "Hit(s) of Interest"

        # Define the treatment and control names
        treatment, control = comparison.split("_vs_")

        # Get the MAGeCK RRA results
        source = MAGeCK_RRA_data(
            comparison=comparison,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )
        # Use numpy.where to conditionally assign the LFC based on the sign of beta
        source["default_rank"] = source["neg|lfc"].rank(ascending=True)

        # Use numpy.where to conditionally assign the FDR based on the sign of LFC
        source["selected_fdr"] = np.where(
            source["neg|lfc"] < 0, source["neg|fdr"], source["pos|fdr"]
        )

        # Set the significant label to non-significant by default
        source["significant"] = non_significant_label

        source = filter_by_threshold(
            source,
            "neg|lfc",
            mageck_rra_lfc,
            mageck_rra_lfc_orientation,
            significant_label,
            "selected_fdr",
            fdr_cutoff,
        )

        # Highlight the elements of interest
        source.loc[source.id.isin(elements), "significant"] = highlight_label

        # Define the domain and range for the color scale
        domain = [significant_label, non_significant_label, highlight_label]
        range_ = [sig, non_sig, "blue"]

        # Define the color scale
        source["significant"] = pd.Categorical(
            source["significant"],
            categories=[non_significant_label, significant_label, highlight_label],
            ordered=True,
        )

        # Order rows by label 'significant' to have the highlighted genes on top
        source = source.sort_values(by="significant", ascending=True)

        # Create the snake plot
        chart = (
            alt.Chart(source, title="MAGeCK RRA (%s)" % comparison)
            .mark_circle(size=60)
            .encode(
                x=alt.X("default_rank:Q", axis=alt.Axis(title="Rank")),
                y=alt.Y(
                    "neg|lfc:Q", axis=alt.Axis(title="%s logFold-change" % treatment)
                ),
                tooltip=[
                    "id",
                    "num",
                    "neg|lfc",
                    "neg|p-value",
                    "neg|fdr",
                    "neg|score",
                    "pos|lfc",
                    "pos|p-value",
                    "pos|fdr",
                    "pos|score",
                    "selected_fdr",
                    "significant",
                    "neg|rank",
                    "default_rank",
                ],
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
            )
            .properties(width=800, height=400)
            .interactive()
        )
        line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule().encode(y="y")
        text = (
            alt.Chart(source.query("significant == 'Hit(s) of Interest'"))
            .mark_text(dy=10, dx=20, color="blue")
            .encode(
                x=alt.X("default_rank:Q"), y=alt.Y("neg|lfc:Q"), text=alt.Text("id")
            )
        )

        chart = (
            (chart + line + text).configure_axis(grid=False).configure_view(stroke=None)
        )
        display(chart)

    def _CRISPhieRmix_snake_plot(
        comparison,
        fdr_cutoff,
        non_sig,
        sig,
        results_directory,
        tools_available,
        elements,
    ):
        # Define the tool
        tool = "CRISPhieRmix"

        # Define the mean log2FoldChange and mean log2FoldChange orientation
        crisphiermix_logfc = float(crisphiermix_logfc_widget.value)
        crisphiermix_logfc_orientation = crisphiermix_logfc_orientation_widget.value

        # Define the significant and non-significant labels
        significant_label = f"FDR < {fdr_cutoff} and logFC {get_pretty_orientation(crisphiermix_logfc_orientation)} {crisphiermix_logfc}"
        non_significant_label = f"FDR ≥ {fdr_cutoff} or logFC {get_pretty_orientation(get_reverse_orientation(crisphiermix_logfc_orientation))} {crisphiermix_logfc}"

        # Define the highlight label
        highlight_label = "Hit(s) of Interest"

        # Define the treatment and control names
        treatment, control = comparison.split("_vs_")

        # Get the CRISPhieRmix results
        source = CRISPhieRmix_data(
            comparison=comparison,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )

        # Compute the default rank of the mean log2FoldChange
        source["default_rank"] = source["mean_log2FoldChange"].rank(method="dense")

        source["significant"] = non_significant_label

        source = filter_by_threshold(
            source,
            "mean_log2FoldChange",
            crisphiermix_logfc,
            crisphiermix_logfc_orientation,
            significant_label,
            "FDR",
            fdr_cutoff,
        )

        # Highlight the elements of interest
        source.loc[source.gene.isin(elements), "significant"] = highlight_label

        # Define the domain and range for the color scale
        domain = [significant_label, non_significant_label, highlight_label]
        range_ = [sig, non_sig, "blue"]

        # Define the color scale
        source["significant"] = pd.Categorical(
            source["significant"],
            categories=[non_significant_label, significant_label, highlight_label],
            ordered=True,
        )

        # Order rows by label 'significant' to have the highlighted genes on top
        source = source.sort_values(by="significant", ascending=True)

        # Create the snake plot
        chart = (
            alt.Chart(source, title=f"CRISPhieRmix ({comparison})")
            .mark_circle(size=60)
            .encode(
                x=alt.X("default_rank:Q", axis=alt.Axis(title="Rank")),
                y=alt.Y(
                    "mean_log2FoldChange:Q",
                    axis=alt.Axis(title=f"{treatment} sgRNAs log2FoldChange average"),
                ),
                tooltip=[
                    "gene",
                    "locfdr",
                    "FDR",
                    "significant",
                    "default_rank",
                    "mean_log2FoldChange",
                ],
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
            )
            .properties(width=800, height=400)
            .interactive()
        )
        line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule().encode(y="y")
        text = (
            alt.Chart(source.query("significant == 'Hit(s) of Interest'"))
            .mark_text(dy=10, dx=20, color="blue")
            .encode(
                x=alt.X("default_rank:Q"),
                y=alt.Y("mean_log2FoldChange:Q"),
                text=alt.Text("gene"),
            )
        )
        chart = (
            (chart + line + text).configure_axis(grid=False).configure_view(stroke=None)
        )
        display(chart)

    def _directional_scoring_method_snake_plot(
        comparison,
        fdr_cutoff,
        non_sig,
        sig,
        results_directory,
        tools_available,
        elements,
    ):
        # Define the tool
        tool = "directional_scoring_method"

        # Define the score and score orientation
        directional_scoring_method_score = float(
            directional_scoring_method_score_widget.value
        )
        directional_scoring_method_score_orientation = (
            directional_scoring_method_score_orientation_widget.value
        )

        # Define the significant and non-significant labels
        significant_label = f"score {get_pretty_orientation(directional_scoring_method_score_orientation)} {directional_scoring_method_score}"
        non_significant_label = f"score {get_pretty_orientation(get_reverse_orientation(directional_scoring_method_score_orientation))} {directional_scoring_method_score}"

        # Define the highlight label
        highlight_label = "Hit(s) of Interest"

        # Define the treatment and control names
        treatment, control = comparison.split("_vs_")

        # Get the Directional Scoring Method results
        source = directional_scoring_method_data(
            comparison=comparison,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )

        # Compute the default rank of the score
        source["default_rank"] = source["score"].rank(method="first")

        # Set the significant label to non-significant by default
        source["significant"] = non_significant_label

        # Filter the data based on the score threshold
        source = filter_by_threshold(
            source,
            "score",
            directional_scoring_method_score,
            directional_scoring_method_score_orientation,
            significant_label,
        )

        # Highlight the elements of interest
        source.loc[source.Gene.isin(elements), "significant"] = highlight_label

        # Define the domain and range for the color scale
        domain = [significant_label, non_significant_label, highlight_label]
        range_ = [sig, non_sig, "blue"]

        # Define the color scale
        source["significant"] = pd.Categorical(
            source["significant"],
            categories=[non_significant_label, significant_label, highlight_label],
            ordered=True,
        )

        # Order rows by label 'significant' to have the highlighted genes on top
        source = source.sort_values(by="significant", ascending=True)

        # Create the snake plot
        chart = (
            alt.Chart(source, title=f"Directional Scoring Method ({comparison})")
            .mark_circle(size=60)
            .encode(
                x=alt.X("default_rank:Q", axis=alt.Axis(title="Rank")),
                y=alt.Y("score:Q", axis=alt.Axis(title=f"{treatment} score")),
                tooltip=[
                    "Gene",
                    "up",
                    "down",
                    "n",
                    "significant",
                    "score",
                    "default_rank",
                ],
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Guides threshold:"),
                ),
            )
            .properties(width=800, height=400)
            .interactive()
        )
        line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule().encode(y="y")
        text = (
            alt.Chart(source.query("significant == 'Hit(s) of Interest'"))
            .mark_text(dy=10, dx=20, color="blue")
            .encode(
                x=alt.X("default_rank:Q"), y=alt.Y("score:Q"), text=alt.Text("Gene")
            )
        )
        chart = (
            (chart + line + text).configure_axis(grid=False).configure_view(stroke=None)
        )
        display(chart)

    def _SSREA_like_snake_plot(
        comparison,
        fdr_cutoff,
        non_sig,
        sig,
        results_directory,
        tools_available,
        elements,
    ):
        # Define the tool
        tool = "SSREA"

        # Define the NES and NES orientation
        ssrea_nes = float(ssrea_nes_widget.value)
        ssrea_nes_orientation = ssrea_nes_orientation_widget.value

        # Define the significant and non-significant labels
        significant_label = f"FDR < {fdr_cutoff} and NES {get_pretty_orientation(ssrea_nes_orientation)} {ssrea_nes}"
        non_significant_label = f"FDR ≥ {fdr_cutoff} or NES {get_pretty_orientation(get_reverse_orientation(ssrea_nes_orientation))} {ssrea_nes}"

        # Define the highlight label
        highlight_label = "Hit(s) of Interest"

        # Define the treatment and control names
        treatment, control = comparison.split("_vs_")

        # Get the SSREA results
        source = SSREA_like_data(
            comparison=comparison,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )

        # Compute the default rank of the NES
        source["default_rank"] = source[["NES"]].rank(method="dense")

        # Set the significant label to non-significant by default
        source["significant"] = non_significant_label

        # Filter the data based on the NES threshold
        source = filter_by_threshold(
            source,
            "NES",
            ssrea_nes,
            ssrea_nes_orientation,
            significant_label,
            "padj",
            fdr_cutoff,
        )

        # Highlight the elements of interest
        source.loc[source.pathway.isin(elements), "significant"] = highlight_label

        # Define the domain and range for the color scale
        domain = [significant_label, non_significant_label, highlight_label]
        range_ = [sig, non_sig, "blue"]

        # Define the color scale
        source["significant"] = pd.Categorical(
            source["significant"],
            categories=[non_significant_label, significant_label, highlight_label],
            ordered=True,
        )

        # Order rows by label 'significant' to have the highlighted genes on top
        source = source.sort_values(by="significant", ascending=True)

        # Create the snake plot
        chart = (
            alt.Chart(source, title=f"SSREA method ({comparison})")
            .mark_circle(size=60)
            .encode(
                x=alt.X("default_rank:Q", axis=alt.Axis(title="Rank")),
                y=alt.Y("NES:Q", axis=alt.Axis(title=f"{treatment} NES")),
                tooltip=[
                    "pathway",
                    "pval",
                    "padj",
                    "ES",
                    "NES",
                    "significant",
                    "size",
                    "default_rank",
                ],
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
            )
            .properties(width=800, height=400)
            .interactive()
        )
        line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule().encode(y="y")
        text = (
            alt.Chart(source.query("significant == 'Hit(s) of Interest'"))
            .mark_text(dy=10, dx=20, color="blue")
            .encode(
                x=alt.X("default_rank:Q"), y=alt.Y("NES:Q"), text=alt.Text("pathway")
            )
        )
        chart = (
            (chart + line + text).configure_axis(grid=False).configure_view(stroke=None)
        )
        display(chart)

    def _BAGEL_snake_plot(
        comparison,
        fdr_cutoff,
        non_sig,
        sig,
        results_directory,
        tools_available,
        elements,
    ):
        # Define the tool
        tool = "BAGEL"

        # Define the BAGEL BF cut-off and orientation
        bagel_bf = float(bagel_bf_widget.value)
        bagel_bf_orientation = bagel_bf_orientation_widget.value

        # Define the significant and non-significant labels
        significant_label = f"FDR < {fdr_cutoff} and BF {get_pretty_orientation(bagel_bf_orientation)} {bagel_bf}"
        non_significant_label = f"FDR ≥ {fdr_cutoff} or BF {get_pretty_orientation(get_reverse_orientation(bagel_bf_orientation))} {bagel_bf}"

        # Define the highlight label
        highlight_label = "Hit(s) of Interest"

        # Define the treatment and control names
        treatment, control = comparison.split("_vs_")

        # Get the BAGEL results
        source = BAGEL_data(
            comparison=comparison,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
            bagel_version=bagel_version,
        )

        # Compute the default rank of the BF scores
        source["default_rank"] = source["BF"].rank(method="dense", ascending=False)

        # Set the significant label to non-significant by default
        source["significant"] = non_significant_label

        # Filter the data based on the BF threshold
        source = filter_by_threshold(
            source,
            "BF",
            bagel_bf,
            bagel_bf_orientation,
            significant_label,
        )

        if bagel_version == 1:
            # rename GENE column to Gene
            source = source.rename(columns={"GENE": "Gene"})

        source.loc[source['Gene'].isin(elements), "significant"] = highlight_label

        # Define the domain and range for the color scale
        domain = [significant_label, non_significant_label, highlight_label]
        range_ = [sig, non_sig, "blue"]

        # Define the color scale
        source["significant"] = pd.Categorical(
            source["significant"],
            categories=[non_significant_label, significant_label, highlight_label],
            ordered=True,
        )

        # Order rows by label 'significant' to have the highlighted genes on top
        source = source.sort_values(by="significant", ascending=True)

        # Create the snake plot
        chart = (
            alt.Chart(source, title=f"BAGEL ({comparison})")
            .mark_circle(size=60)
            .encode(
                x=alt.X("default_rank:Q", axis=alt.Axis(title="Rank")),
                y=alt.Y("BF:Q", axis=alt.Axis(title="Bayesian Factor")),
                tooltip=source.columns.tolist() + ["default_rank"],
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
            )
            .properties(width=800, height=400)
            .interactive()
        )
        line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule().encode(y="y")
        text = (
            alt.Chart(source.query("significant == 'Hit(s) of Interest'"))
            .mark_text(dy=10, dx=20, color="blue")
            .encode(x=alt.X("default_rank:Q"), y=alt.Y("BF:Q"), text=alt.Text("Gene"))
        )
        chart = (
            (chart + line + text).configure_axis(grid=False).configure_view(stroke=None)
        )
        display(chart)

    def _plot(event):
        tool = tools_widget.value
        comparison = comparisons_widget.value
        fdr_cutoff = fdr_widget.value
        non_sig = color_non_widget.value
        sig = color_sig_widget.value
        elements = element.value.split(",")
        treatment = comparison.split("_vs_")[0]
        control = comparison.split("_vs_")[1]
        at_least_one_tool = False
        if "MAGeCK_RRA" in tool:
            at_least_one_tool = True
            _MAGeCK_RRA_snake_plot(
                comparison,
                fdr_cutoff,
                non_sig,
                sig,
                results_directory,
                tools_available,
                elements,
            )
            download(
                tools_available, tool="MAGeCK_RRA", treatment=treatment, control=control
            )
        if "MAGeCK_MLE" in tool:
            at_least_one_tool = True
            _MAGeCK_MLE_snake_plot(
                comparison,
                fdr_cutoff,
                non_sig,
                sig,
                results_directory,
                tools_available,
                elements,
            )
            download(
                tools_available, tool="MAGeCK_MLE", treatment=treatment, control=control
            )
        if "CRISPhieRmix" in tool:
            at_least_one_tool = True
            _CRISPhieRmix_snake_plot(
                comparison,
                fdr_cutoff,
                non_sig,
                sig,
                results_directory,
                tools_available,
                elements,
            )
            download(
                tools_available,
                tool="CRISPhieRmix",
                treatment=treatment,
                control=control,
            )
        if "directional_scoring_method" in tool:
            at_least_one_tool = True
            _directional_scoring_method_snake_plot(
                comparison,
                fdr_cutoff,
                non_sig,
                sig,
                results_directory,
                tools_available,
                elements,
            )
            download(
                tools_available,
                tool="directional_scoring_method",
                treatment=treatment,
                control=control,
            )
        if "SSREA" in tool:
            at_least_one_tool = True
            _SSREA_like_snake_plot(
                comparison,
                fdr_cutoff,
                non_sig,
                sig,
                results_directory,
                tools_available,
                elements,
            )
            download(
                tools_available, tool="SSREA", treatment=treatment, control=control
            )
        if "BAGEL" in tool:
            at_least_one_tool = True
            _BAGEL_snake_plot(
                comparison,
                fdr_cutoff,
                non_sig,
                sig,
                results_directory,
                tools_available,
                elements,
            )
            download(
                tools_available, tool="BAGEL", treatment=treatment, control=control
            )
        if not at_least_one_tool:
            # print("Choose a tool.")
            display_warning("Choose at least one tool.")

    display(tools_widget)
    display(element)
    display(comparisons_widget)
    # Display a text widget to delimite the filter for each tool
    display(HTML("<h3>Highlight features:</h3>"))
    display(fdr_widget)

    if "BAGEL" in tools:
        display(HTML("<h4>BAGEL:</h4>"))
        display(widgets.HBox([bagel_bf_widget, bagel_bf_orientation_widget]))

    if "SSREA" in tools:
        display(HTML("<h4>SSREA:</h4>"))
        display(widgets.HBox([ssrea_nes_widget, ssrea_nes_orientation_widget]))

    if "directional_scoring_method" in tools:
        display(HTML("<h4>DSM:</h4>"))
        display(
            widgets.HBox(
                [
                    directional_scoring_method_score_widget,
                    directional_scoring_method_score_orientation_widget,
                ]
            )
        )

    if "CRISPhieRmix" in tools:
        display(HTML("<h4>CRISPhieRmix:</h4>"))
        display(
            widgets.HBox(
                [crisphiermix_logfc_widget, crisphiermix_logfc_orientation_widget]
            )
        )

    if "MAGeCK_MLE" in tools:
        display(HTML("<h4>MAGeCK MLE:</h4>"))
        display(
            widgets.HBox([mageck_mle_beta_widget, mageck_mle_beta_orientation_widget])
        )

    if "MAGeCK_RRA" in tools:
        display(HTML("<h4>MAGeCK RRA:</h4>"))
        display(
            widgets.HBox([mageck_rra_lfc_widget, mageck_rra_lfc_orientation_widget])
        )

    display(HTML("<h3>Color:</h3>"))
    display(color_sig_widget)
    display(color_non_widget)

    button = widgets.Button(description="Click!")

    display(button)

    button.on_click(_plot)


def show_sgRNA_counts(token):
    config = "./config/%s.yaml" % token
    content = open_yaml(config)
    cts_file = "results/%s/normalized.filtered.counts.txt" % token
    cts = pd.read_csv(cts_file, sep="\t")
    cts_columns = [col for col in cts.columns.tolist() if not col in ["sgRNA", "Gene"]]
    cts = pd.melt(cts, id_vars=["sgRNA", "Gene"])
    genes_list = list(set(cts.Gene.tolist()))

    color_schemes = [
        "blueorange",
        "brownbluegreen",
        "purplegreen",
        "pinkyellowgreen",
        "purpleorange",
        "redblue",
        "redgrey",
        "redyellowblue",
        "redyellowgreen",
        "spectral",
    ]

    element = widgets.Combobox(
        placeholder="Choose one",
        options=genes_list,
        description="Element:",
        value=genes_list[0],
        ensure_option=False,
    )
    conditions = widgets.TagsInput(
        value=cts_columns, allowed_tags=cts_columns, allow_duplicates=False
    )
    color_scheme_widget = widgets.Dropdown(
        placeholder="Choose one",
        options=color_schemes,
        description="Color scheme:",
        value=color_schemes[0],
    )
    reverse_widget = widgets.Checkbox(value=False, description="Reverse color scheme?")

    display(widgets.VBox([element, conditions, color_scheme_widget, reverse_widget]))

    button = widgets.Button(description="Show!")
    display(button)

    def on_button_clicked(b):
        genes = element.value.split(",")
        for gene in genes:
            if not gene in list(cts.Gene):
                # print("%s was not found in %s." % (gene, cts_file))
                display_warning(f"Gene '{gene}' was not found in {cts_file}.")
            else:
                gene_cts = cts.loc[cts.Gene == gene]
                gene_cts = gene_cts.loc[gene_cts.variable.isin(conditions.value)]
                z_scores = gene_cts.groupby(["sgRNA"]).value.transform(
                    lambda x: zscore(x, ddof=1)
                )
                gene_cts["z-score"] = z_scores
                chart = (
                    alt.Chart(
                        gene_cts, title="%s normalized read counts heatmap" % gene
                    )
                    .mark_rect()
                    .encode(
                        x=alt.X(
                            "variable",
                            axis=alt.Axis(title="Replicate"),
                            sort=conditions.value,
                        ),
                        y=alt.Y("sgRNA", axis=alt.Axis(title="sgRNA")),
                        color=alt.Color(
                            "z-score:Q",
                            scale=alt.Scale(
                                scheme=color_scheme_widget.value,
                                reverse=reverse_widget.value,
                            ),
                        ),
                        tooltip=["sgRNA", "variable", "value", "z-score"],
                    )
                    .interactive()
                )
                display(chart)

    button.on_click(on_button_clicked)


def show_sgRNA_counts_lines(token):
    config = "config/%s.yaml" % token
    content = open_yaml(config)
    cts_file = "results/%s/normalized.filtered.counts.txt" % token
    cts = pd.read_csv(cts_file, sep="\t")
    design_file = content["tsv_file"]
    design = pd.read_csv(design_file, sep="\t")

    conditions_list_init = list(design.condition)
    conditions_list = []
    for condition in conditions_list_init:
        if not condition in conditions_list:
            conditions_list.append(condition)

    genes_list = list(set(pd.melt(cts, id_vars=["sgRNA", "Gene"]).Gene.tolist()))
    element = widgets.Combobox(
        placeholder="Choose one",
        options=genes_list,
        description="Element:",
        value=genes_list[0],
        ensure_option=False,
    )
    conditions = widgets.TagsInput(
        value=conditions_list, allowed_tags=conditions_list, allow_duplicates=False
    )
    button = widgets.Button(description="Show!")

    display(element)
    display(conditions)
    display(button)

    def show_plot(source, sort_cols, gene):
        selection = alt.selection_multi(fields=["sgRNA"], bind="legend")

        line = (
            alt.Chart(source)
            .mark_line()
            .encode(
                x=alt.X("condition:O", sort=sort_cols),
                y="mean_value:Q",
                color=alt.Color("sgRNA:N"),
                opacity=alt.condition(selection, alt.value(1), alt.value(0.0)),
                tooltip=["sgRNA:N", "mean_value:Q"],
            )
            .transform_aggregate(
                mean_value="mean(value)", groupby=["condition", "sgRNA"]
            )
            .add_selection(selection)
            .properties(width=600)
        )

        band = (
            alt.Chart(source)
            .mark_errorband(extent="stdev")
            .encode(
                x=alt.X("condition:O", sort=sort_cols),
                y="value:Q",
                color=alt.Color("sgRNA:N"),
                opacity=alt.condition(selection, alt.value(0.5), alt.value(0.0)),
            )
        )

        chart = line + band
        display(chart.interactive())

    def on_button_clicked(b):
        # Read normalized counts file
        cts = pd.read_csv(cts_file, sep="\t")
        # Transform dataframe: one value per line
        cts = pd.melt(cts, id_vars=["sgRNA", "Gene"])
        genes = element.value.split(",")
        for gene in genes:
            if not gene in list(cts.Gene):
                display_warning(f"Gene '{gene}' was not found in {cts_file}.")
            else:
                sort_cols = conditions.value
                gene_cts = cts.loc[cts.Gene == gene]
                source = gene_cts
                source = pd.merge(
                    source, design, left_on="variable", right_on="replicate"
                )
                boolean_series = source.condition.isin(sort_cols)
                source = source[boolean_series]
                # source["replicate_group"] = source["replicate"].str.extract(
                #     r"([1-9]+)$"
                # )
                # source = pd.pivot_table(
                #     source,
                #     values="value",
                #     index=["sgRNA", "replicate"],
                #     columns=["condition"],
                # )
                source = source.reset_index()
                show_plot(source, sort_cols, gene)

    button.on_click(on_button_clicked)


def regions2genes(token, se):
    annotation_file = f"resources/{token}/annotation_ROSE_REGION_TO_GENE.txt"
    # print(os.path.exists(annotation_file), annotation_file)
    if os.path.exists(annotation_file):
        annotation_table = pd.read_table(
            annotation_file, sep="\t", skiprows=1, header=None
        )
        annotation_table.columns = [
            "Name",
            "Chromosome",
            "Start",
            "End",
            "OVERLAP_GeneS",
            "PROXIMAL_GeneS",
            "CLOSEST_Gene",
            "L",
            "none2",
        ]
        # print(annotation_table.head())
        annotation_table = annotation_table.loc[annotation_table.Name.isin(se)]
        annotation_table["OVERLAP_GeneS"] = annotation_table["OVERLAP_GeneS"].fillna(0)
        annotation_table["PROXIMAL_GeneS"] = annotation_table["PROXIMAL_GeneS"].fillna(
            0
        )
        annotation_table["CLOSEST_Gene"] = annotation_table["CLOSEST_Gene"].fillna(0)
        overlap = [x for x in annotation_table["OVERLAP_GeneS"].values if x != 0]
        proximal = [x for x in annotation_table["PROXIMAL_GeneS"].values if x != 0]
        closest = [x for x in annotation_table["CLOSEST_Gene"].values if x != 0]
        s2g = []
        for genes_set in [overlap, proximal, closest]:
            for i in range(len(genes_set)):
                genes = genes_set[i].split(",")
                for gene in genes:
                    if not gene in s2g:
                        s2g.append(gene)
    else:
        s2g = se
    # print(s2g)
    return s2g


def tool_results_by_element(results_directory, tools_available, token, bagel_version):
    def get_controls(results_directory, tools_available, tool):
        comparisons_list = os.listdir(os.path.join(results_directory, tool))
        ctrs = list(
            set(map(lambda x: x.split("_vs_")[1], list(tools_available[tool].keys())))
        )
        return ctrs

    def get_tool_results(results_directory, tools_available, tool):
        if tool == "CRISPhieRmix":
            result = CRISPhieRmix_data(
                comparison="",
                control=control.value,
                tool=tool,
                results_directory=results_directory,
                tools_available=tools_available,
            )
        elif tool == "MAGeCK_MLE":
            result = MAGeCK_MLE_data(
                comparison="",
                control=control.value,
                tool=tool,
                results_directory=results_directory,
                tools_available=tools_available,
            )
        elif tool == "MAGeCK_RRA":
            result = MAGeCK_RRA_data(
                comparison="",
                control=control.value,
                tool=tool,
                results_directory=results_directory,
                tools_available=tools_available,
            )
        elif tool == "SSREA":
            result = SSREA_like_data(
                comparison="",
                control=control.value,
                tool=tool,
                results_directory=results_directory,
                tools_available=tools_available,
            )
        elif tool == "directional_scoring_method":
            result = directional_scoring_method_data(
                comparison="",
                control=control.value,
                tool=tool,
                results_directory=results_directory,
                tools_available=tools_available,
            )
        if tool == "BAGEL":
            result = BAGEL_data(
                comparison="",
                control=control.value,
                tool=tool,
                results_directory=results_directory,
                tools_available=tools_available,
                bagel_version=bagel_version
            )
        return result

    def get_genes_list(results_directory, tools_available, tool):
        result = get_tool_results(results_directory, tools_available, tool)
        if tool == "CRISPhieRmix":
            elements_list = list(set(result.gene))
        elif tool == "MAGeCK_MLE":
            elements_list = list(set(result.Gene))
        elif tool == "MAGeCK_RRA":
            elements_list = list(set(result.id))
        elif tool == "BAGEL":
            elements_list = list(set(result.Gene))
        elif tool == "SSREA":
            elements_list = list(set(result.pathway))
        elif tool == "directional_scoring_method":
            elements_list = list(set(result.Gene))
        return elements_list

    def update_genes_list(update):
        tools = tools_available.keys()
        for tool in tools:
            if tool == "CRISPhieRmix":
                result = CRISPhieRmix_data(
                    comparison="",
                    control=control.value,
                    tool=tool,
                    results_directory=results_directory,
                    tools_available=tools_available,
                )
                gene_var = "gene"
                break
            if tool == "MAGeCK_MLE":
                result = MAGeCK_MLE_data(
                    comparison="",
                    control=control.value,
                    tool=tool,
                    results_directory=results_directory,
                    tools_available=tools_available,
                )
                gene_var = "Gene"
                break
            if tool == "MAGeCK_RRA":
                result = MAGeCK_RRA_data(
                    comparison="",
                    control=control.value,
                    tool=tool,
                    results_directory=results_directory,
                    tools_available=tools_available,
                )
                gene_var = "ig"
                break
            if tool == "BAGEL":
                result = BAGEL_data(
                    comparison="",
                    control=control.value,
                    tool=tool,
                    results_directory=results_directory,
                    tools_available=tools_available,
                    bagel_version=bagel_version
                )
                gene_var = "Gene"
                break
            if tool == "SSREA":
                result = SSREA_like_data(
                    comparison="",
                    control=control.value,
                    tool=tool,
                    results_directory=results_directory,
                    tools_available=tools_available,
                )
                gene_var = "pathway"
                break
            if tool == "directional_scoring_method":
                result = directional_scoring_method_data(
                    comparison="",
                    control=control.value,
                    tool=tool,
                    results_directory=results_directory,
                    tools_available=tools_available,
                )
                gene_var = "Gene"
                break

        elements_list = result[gene_var].to_list()
        gene.options = list(set(elements_list))
        gene.value = elements_list[0]

    # Set configuration
    config = f"config/{token}.yaml"

    # Open configuration file
    content = open_yaml(config)
    design_file = content["tsv_file"]
    design = pd.read_csv(design_file, sep="\t")

    conditions_list_init = list(design.condition)
    conditions_list = []
    for condition in conditions_list_init:
        if not condition in conditions_list:
            conditions_list.append(condition)

    tools_list = [tool for tool in list(tools_available.keys()) if tool != "DESeq2"]

    ctrs = get_controls(results_directory, tools_available, tools_list[0])
    control = widgets.Dropdown(
        options=ctrs, value=ctrs[0], description="Control:", disabled=False
    )
    control.observe(update_genes_list, "value")

    conditions = widgets.TagsInput(
        value=conditions_list, allowed_tags=conditions_list, allow_duplicates=False
    )

    fdr_widget = widgets.FloatSlider(
        min=0.0, max=1.0, step=0.01, value=0.05, description="FDR cut-off:"
    )

    # Add a widget to define the BAGEL minimum BF cut-off. Default is 0. No minimum and maximum values are defined. Using FloatText widget.
    bagel_bf_widget = widgets.FloatText(
        value=0,
        description="BAGEL BF cut-off:",
        # If BAGEL is not available, the widget is disabled
        disabled="BAGEL" not in tools_list,
        style=style,
    )

    # Add a widget to define the orientation of the BAGEL BF cut-off. Default is ">=".
    bagel_bf_orientation_widget = widgets.Dropdown(
        options=[">=", "<=", "abs() >=", "abs() <="],
        value=">=",
        description="",
        # If BAGEL is not available, the widget is disabled
        disabled="BAGEL" not in tools_list,
        # style=style,
        layout=widgets.Layout(width="75px"),
    )

    # Add a widget to define the SSREA minimum NES cut-off. Default is 0. No minimum and maximum values are defined. Using FloatText widget.
    ssrea_nes_widget = widgets.FloatText(
        value=0,
        description="SSREA NES cut-off:",
        # If SSREA is not available, the widget is disabled
        disabled="SSREA" not in tools_list,
        style=style,
    )

    # Add a widget to define the orientation of the SSREA NES cut-off. Default is ">=". Using Dropdown widget.
    ssrea_nes_orientation_widget = widgets.Dropdown(
        options=[">=", "<=", "abs() >=", "abs() <="],
        value="abs() >=",
        description="",
        # If SSREA is not available, the widget is disabled
        disabled="SSREA" not in tools_list,
        # style=style,
        layout=widgets.Layout(width="75px"),
    )

    # Add a widget to define the directional_scoring_method minimum score cut-off. Default is 0. No minimum and maximum values are defined. Using FloatText widget.
    directional_scoring_method_score_widget = widgets.FloatText(
        value=0,
        description="DSM score cut-off:",
        # If directional_scoring_method is not available, the widget is disabled
        disabled="directional_scoring_method" not in tools_list,
        style=style,
    )

    # Add a widget to define the orientation of the directional_scoring_method score cut-off. Default is ">=". Using Dropdown widget.
    directional_scoring_method_score_orientation_widget = widgets.Dropdown(
        options=[">=", "<=", "abs() >=", "abs() <=", "abs() >", "abs() <"],
        value="abs() >",
        description="",
        # If directional_scoring_method is not available, the widget is disabled
        disabled="directional_scoring_method" not in tools_list,
        # style=style,
        layout=widgets.Layout(width="75px"),
    )

    # Add a widget to define the CRISPhieRmix minimum mean logFC cut-off. Default is 0. No minimum and maximum values are defined. Using FloatText widget.
    crisphiermix_logfc_widget = widgets.FloatText(
        value=0,
        description="CRISPhieRmix logFC cut-off:",
        # If CRISPhieRmix is not available, the widget is disabled
        disabled="CRISPhieRmix" not in tools_list,
        style=style,
    )

    # Add a widget to define the orientation of the CRISPhieRmix mean logFC cut-off. Default is ">=". Using Dropdown widget.
    crisphiermix_logfc_orientation_widget = widgets.Dropdown(
        options=[">=", "<=", "abs() >=", "abs() <="],
        value="abs() >=",
        description="",
        # If CRISPhieRmix is not available, the widget is disabled
        disabled="CRISPhieRmix" not in tools_list,
        # style=style,
        layout=widgets.Layout(width="75px"),
    )

    # Add a widget to define the MAGeCK_MLE minimum beta cut-off. Default is 0. No minimum and maximum values are defined. Using FloatText widget.
    mageck_mle_beta_widget = widgets.FloatText(
        value=0,
        description="MAGeCK MLE beta cut-off:",
        # If MAGeCK_MLE is not available, the widget is disabled
        disabled="MAGeCK_MLE" not in tools_list,
        style=style,
    )

    # Add a widget to define the orientation of the MAGeCK_MLE beta cut-off. Default is ">=". Using Dropdown widget.
    mageck_mle_beta_orientation_widget = widgets.Dropdown(
        options=[">=", "<=", "abs() >=", "abs() <="],
        value="abs() >=",
        description="",
        # If MAGeCK_MLE is not available, the widget is disabled
        disabled="MAGeCK_MLE" not in tools_list,
        # style=style,
        layout=widgets.Layout(width="75px"),
    )

    # Add a widget to define the MAGeCK_RRA minimum LFC cut-off. Default is 0. No minimum and maximum values are defined. Using FloatText widget.
    mageck_rra_lfc_widget = widgets.FloatText(
        value=0,
        description="MAGeCK RRA LFC cut-off:",
        # If MAGeCK_RRA is not available, the widget is disabled
        disabled="MAGeCK_RRA" not in tools_list,
        style=style,
    )

    # Add a widget to define the orientation of the MAGeCK_RRA LFC cut-off. Default is ">=". Using Dropdown widget.
    mageck_rra_lfc_orientation_widget = widgets.Dropdown(
        options=[">=", "<=", "abs() >=", "abs() <="],
        value="abs() >=",
        description="",
        # If MAGeCK_RRA is not available, the widget is disabled
        disabled="MAGeCK_RRA" not in tools_list,
        # style=style,
        layout=widgets.Layout(width="75px"),
    )

    elements_list = get_genes_list(results_directory, tools_available, tools_list[0])
    gene = widgets.Combobox(
        placeholder="Choose one",
        options=elements_list,
        description="Element:",
        value=elements_list[0],
        ensure_option=False,
    )

    # Display all widgets
    display(control)
    display(gene)
    display(fdr_widget)

    if "BAGEL" in tools_list:
        display(widgets.HBox([bagel_bf_widget, bagel_bf_orientation_widget]))

    if "SSREA" in tools_list:
        display(widgets.HBox([ssrea_nes_widget, ssrea_nes_orientation_widget]))

    if "directional_scoring_method" in tools_list:
        display(
            widgets.HBox(
                [
                    directional_scoring_method_score_widget,
                    directional_scoring_method_score_orientation_widget,
                ]
            )
        )

    if "CRISPhieRmix" in tools_list:
        display(
            widgets.HBox(
                [crisphiermix_logfc_widget, crisphiermix_logfc_orientation_widget]
            )
        )

    if "MAGeCK_MLE" in tools_list:
        display(
            widgets.HBox([mageck_mle_beta_widget, mageck_mle_beta_orientation_widget])
        )

    if "MAGeCK_RRA" in tools_list:
        display(
            widgets.HBox([mageck_rra_lfc_widget, mageck_rra_lfc_orientation_widget])
        )

    def CRISPhieRmix_results(result, fdr_cutoff, control, gene, sort_cols):
        # Define the orientation of the logFC cut-off and its value
        crisphiermix_logfc_orientation = crisphiermix_logfc_orientation_widget.value
        crisphiermix_logfc = crisphiermix_logfc_widget.value

        # Define the significant and non-significant labels
        significant_label = f"FDR < {fdr_cutoff} and logFC {get_pretty_orientation(crisphiermix_logfc_orientation)} {crisphiermix_logfc}"
        non_significant_label = f"FDR ≥ {fdr_cutoff} or logFC {get_pretty_orientation(get_reverse_orientation(crisphiermix_logfc_orientation))} {crisphiermix_logfc}"

        #  Define default values for the significant column
        result["significant"] = non_significant_label

        # Filter the data based on the logFC threshold
        result = filter_by_threshold(
            result,
            "mean_log2FoldChange",
            crisphiermix_logfc,
            crisphiermix_logfc_orientation,
            significant_label,
        )

        # Add a new row to the dataframe to display the baseline
        new_row = {
            "gene": gene,
            "condition": control,
            "significant": "Baseline",
            "locfdr": 1,
            "mean_log2FoldChange": 0,
        }

        # Append the new row to the dataframe
        result = result.append(new_row, ignore_index=True)

        # Filter the dataframe to keep only the conditions in sort_cols
        res = result.loc[result.gene == gene]

        # filter res to keep 'condition' in sort_cols
        if control not in sort_cols:
            sort_cols.append(control)

        # Filter the dataframe to keep only the conditions in sort_cols
        res = res[res["condition"].isin(sort_cols)]

        # Define the domain and range for the color scale
        domain = [significant_label, non_significant_label, "Baseline"]
        range_ = ["red", "grey", "black"]

        # Create the plot
        plot = (
            alt.Chart(res)
            .transform_fold(sort_cols)
            .mark_circle(size=60)
            .mark_point(filled=True, size=100, opacity=1.0)
            .encode(
                y=alt.Y(
                    "mean_log2FoldChange",
                    axis=alt.Axis(title="sgRNAs log2FoldChange average"),
                ),
                x=alt.X("condition:N", sort=sort_cols),
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                tooltip=["gene", "locfdr", "significant", "mean_log2FoldChange"],
            )
            .properties(title=gene + " (CRISPhieRmix)", width=100)
        )
        return plot

    def MAGeCK_RRA_results(result, fdr_cutoff, control, gene, sort_cols):
        # Define the threshold for the score and its orientation
        mageck_rra_lfc_orientation = mageck_rra_lfc_orientation_widget.value
        mageck_rra_lfc = mageck_rra_lfc_widget.value

        # Define the significant and non-significant labels
        significant_label = f"FDR < {fdr_cutoff} and LFC {get_pretty_orientation(mageck_rra_lfc_orientation)} {mageck_rra_lfc}"
        non_significant_label = f"FDR ≥ {fdr_cutoff} or LFC {get_pretty_orientation(get_reverse_orientation(mageck_rra_lfc_orientation))} {mageck_rra_lfc}"

        # Use numpy.where to conditionally assign the FDR based on the sign of LFC
        result["selected_fdr"] = np.where(
            result["neg|lfc"] < 0, result["neg|fdr"], result["pos|fdr"]
        )

        # Define default values for the significant column
        result["significant"] = non_significant_label

        # Filter the data based on the LFC threshold
        result = filter_by_threshold(
            result,
            "neg|lfc",
            mageck_rra_lfc,
            mageck_rra_lfc_orientation,
            significant_label,
        )

        new_row = {
            "id": gene,
            "condition": control,
            "significant": "Baseline",
            "neg|fdr": 1,
            "neg|lfc": 0,
            "neg|p-value": 1,
            "pos|fdr": 1,
            "pos|lfc": 0,
            "neg|score": 1,
            "pos|score": 1,
            "pos|p-value": 1,
        }

        result = result.append(new_row, ignore_index=True)

        res = result.loc[result.id == gene]

        res["score"] = np.where(res["neg|lfc"] < 0, res["neg|score"], res["pos|score"])

        # filter res to keep 'condition' in sort_cols
        if control not in sort_cols:
            sort_cols.append(control)
        res = res[res["condition"].isin(sort_cols)]
        domain = [significant_label, non_significant_label, "Baseline"]
        range_ = ["red", "grey", "black"]
        plot = (
            alt.Chart(res)
            .mark_circle(size=60)
            .mark_point(filled=True, size=100, opacity=1.0)
            .encode(
                y=alt.Y(f"score", axis=alt.Axis(title="RRA Score")),
                x=alt.X("condition:N", sort=sort_cols),
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                tooltip=[
                    "id",
                    "neg|score",
                    "neg|lfc",
                    "neg|fdr",
                    "neg|p-value",
                    "pos|lfc",
                    "pos|fdr",
                    "pos|score",
                    "pos|p-value",
                    "score",
                    "significant",
                    "condition",
                ],
            )
            .properties(title=gene + " (MAGeCK RRA)", width=100)
        )
        return plot

    def MAGeCK_MLE_results(result, fdr_cutoff, control, gene, sort_cols):

        # Define the threshold for the beta and its orientation
        mageck_mle_beta_orientation = mageck_mle_beta_orientation_widget.value
        mageck_mle_beta = mageck_mle_beta_widget.value

        # Define the significant and non-significant labels
        significant_label = f"FDR < {fdr_cutoff} and beta {get_pretty_orientation(mageck_mle_beta_orientation)} {mageck_mle_beta}"
        non_significant_label = f"FDR ≥ {fdr_cutoff} or beta {get_pretty_orientation(get_reverse_orientation(mageck_mle_beta_orientation))} {mageck_mle_beta}"

        # Define default values for the significant column
        result["significant"] = non_significant_label

        # Filter the data based on the beta threshold
        result = filter_by_threshold(
            result,
            "beta",
            mageck_mle_beta,
            mageck_mle_beta_orientation,
            significant_label,
        )

        # Add a new row to the dataframe to display the baseline
        new_row = {
            "Gene": gene,
            "condition": control,
            "significant": "Baseline",
            "fdr": 1,
            "beta": 0,
        }
        result = result.append(new_row, ignore_index=True)

        # Filter the dataframe to keep only the conditions in sort_cols
        res = result.loc[result.Gene == gene]

        # Filter the dataframe to keep only the conditions in sort_cols
        if control not in sort_cols:
            sort_cols.append(control)

        res = res[res["condition"].isin(sort_cols)]
        domain = [significant_label, non_significant_label, "Baseline"]
        range_ = ["red", "grey", "black"]
        plot = (
            alt.Chart(res)
            .mark_circle(size=60)
            .mark_point(filled=True, size=100, opacity=1.0)
            .encode(
                y=alt.Y("beta", axis=alt.Axis(title="Beta")),
                x=alt.X("condition:N", sort=sort_cols),
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                tooltip=["Gene", "beta", "significant", "fdr", "condition"],
            )
            .properties(title=gene + " (MAGeCK MLE)", width=100)
        )
        return plot

    def SSREA_like_results(result, fdr_cutoff, control, gene, sort_cols):

        # Define the threshold for the NES and its orientation
        ssrea_nes_orientation = ssrea_nes_orientation_widget.value
        ssrea_nes = ssrea_nes_widget.value

        # Define the significant and non-significant labels
        significant_label = f"FDR < {fdr_cutoff} and NES {get_pretty_orientation(ssrea_nes_orientation)} {ssrea_nes}"
        non_significant_label = f"FDR ≥ {fdr_cutoff} or NES {get_pretty_orientation(get_reverse_orientation(ssrea_nes_orientation))} {ssrea_nes}"

        # Define default values for the significant column
        result["significant"] = non_significant_label

        # Filter the data based on the NES threshold
        result = filter_by_threshold(
            result, "NES", ssrea_nes, ssrea_nes_orientation, significant_label
        )

        # Add a new row to the dataframe to display the baseline
        new_row = {
            "pathway": gene,
            "condition": control,
            "significant": "Baseline",
            "pval": 1,
            "padj": 1,
            "ES": 0,
            "NES": 0,
            "size": None,
        }

        # Append the new row to the dataframe
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.pathway == gene]
        # filter res to keep 'condition' in sort_cols
        if control not in sort_cols:
            sort_cols.append(control)
        res = res[res["condition"].isin(sort_cols)]
        domain = [significant_label, non_significant_label, "Baseline"]
        range_ = ["red", "grey", "black"]
        plot = (
            alt.Chart(res)
            .mark_circle(size=60)
            .mark_point(filled=True, size=100, opacity=1.0)
            .encode(
                y=alt.Y("NES", axis=alt.Axis(title="NES")),
                x=alt.X("condition:N", sort=sort_cols),
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                tooltip=["pathway", "condition", "significant", "pval", "padj", "NES"],
            )
            .properties(title=gene + " (SSREA method)", width=100)
        )
        return plot

    def BAGEL_results(result, fdr_cutoff, control, gene, sort_cols):

        # Define the threshold for the BF and its orientation
        bagel_bf_orientation = bagel_bf_orientation_widget.value
        bagel_bf = bagel_bf_widget.value

        # Define the significant and non-significant labels
        significant_label = (
            f"BF {get_pretty_orientation(bagel_bf_orientation)} {bagel_bf}"
        )
        non_significant_label = f"BF {get_pretty_orientation(get_reverse_orientation(bagel_bf_orientation))} {bagel_bf}"

        # Define default values for the significant column
        result["significant"] = non_significant_label

        # Filter the data based on the BF threshold
        result = filter_by_threshold(
            result, "BF", bagel_bf, bagel_bf_orientation, significant_label
        )

        # rename GENE to Gene
        result = result.rename(columns={"GENE": "Gene"})

        new_row = {
            "Gene": gene,
            "condition": control,
            "significant": "Baseline",
            "BF": 0,
            "FDR": 1,
        }
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.Gene == gene]
        # filter res to keep 'condition' in sort_cols
        if control not in sort_cols:
            sort_cols.append(control)
        res = res[res["condition"].isin(sort_cols)]
        domain = [significant_label, non_significant_label, "Baseline"]
        range_ = ["red", "grey", "black"]
        plot = (
            alt.Chart(res)
            .mark_circle(size=60)
            .mark_point(filled=True, size=100, opacity=1.0)
            .encode(
                y=alt.Y("BF", axis=alt.Axis(title="Bayesian factor")),
                x=alt.X("condition:N", sort=sort_cols),
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                tooltip=res.columns.tolist(),
            )
            .properties(title=gene + " (BAGEL)", width=100)
        )
        return plot

    def directional_scoring_method_results(
        result, fdr_cutoff, control, gene, sort_cols
    ):

        # Define the threshold for the score and its orientation
        directional_scoring_method_score_orientation = (
            directional_scoring_method_score_orientation_widget.value
        )
        directional_scoring_method_score = directional_scoring_method_score_widget.value

        # Define the significant and non-significant labels
        significant_label = f"score {get_pretty_orientation(directional_scoring_method_score_orientation)} {directional_scoring_method_score}"
        non_significant_label = f"score {get_pretty_orientation(get_reverse_orientation(directional_scoring_method_score_orientation))} {directional_scoring_method_score}"

        # Define default values for the significant column
        result["significant"] = non_significant_label

        # Filter the data based on the score threshold
        result = filter_by_threshold(
            result,
            "score",
            directional_scoring_method_score,
            directional_scoring_method_score_orientation,
            significant_label,
        )

        new_row = {
            "Gene": gene,
            "condition": control,
            "significant": "Baseline",
            "down": 0,
            "score": 0,
            "up": 0,
            "n": 0,
            "category": "Baseline",
        }
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.Gene == gene]
        # filter res to keep 'condition' in sort_cols
        if control not in sort_cols:
            sort_cols.append(control)
        res = res[res["condition"].isin(sort_cols)]
        domain = [significant_label, non_significant_label, "Baseline"]
        range_ = ["red", "grey", "black"]
        plot = (
            alt.Chart(res)
            .mark_circle(size=60)
            .mark_point(filled=True, size=100, opacity=1.0)
            .encode(
                y="score",
                x=alt.X("condition:N", sort=sort_cols),
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                tooltip=[
                    "Gene",
                    "condition",
                    "down",
                    "up",
                    "n",
                    "score",
                    "significant",
                ],
            )
            .properties(title=gene + " (Directional Scoring Method)", width=100)
        )
        return plot

    def on_button_clicked(b):
        genes = gene.value.split(",")
        for element in genes:
            chart = alt.hconcat()
            for tool in tools_list:
                result = get_tool_results(results_directory, tools_available, tool)
                if tool == "CRISPhieRmix":
                    plot = CRISPhieRmix_results(
                        result,
                        fdr_widget.value,
                        control.value,
                        element,
                        conditions.value,
                    )
                elif tool == "MAGeCK_MLE":
                    plot = MAGeCK_MLE_results(
                        result,
                        fdr_widget.value,
                        control.value,
                        element,
                        conditions.value,
                    )
                elif tool == "SSREA":
                    plot = SSREA_like_results(
                        result,
                        fdr_widget.value,
                        control.value,
                        element,
                        conditions.value,
                    )
                elif tool == "directional_scoring_method":
                    plot = directional_scoring_method_results(
                        result,
                        fdr_widget.value,
                        control.value,
                        element,
                        conditions.value,
                    )
                elif tool == "MAGeCK_RRA":
                    plot = MAGeCK_RRA_results(
                        result,
                        fdr_widget.value,
                        control.value,
                        element,
                        conditions.value,
                    )
                elif tool == "BAGEL":
                    plot = BAGEL_results(
                        result,
                        fdr_widget.value,
                        control.value,
                        element,
                        conditions.value,
                    )
                chart |= plot
            chart = chart.resolve_legend(
                color="independent", size="independent"
            ).resolve_scale(color="independent", size="independent")
            display(chart)

    button = widgets.Button(description="Show plot")
    display(button)
    button.on_click(on_button_clicked)


def get_enrichr_bases():
    """Get the list of available libraries in EnrichR"""
    link = "https://maayanlab.cloud/Enrichr/datasetStatistics"

    try:
        response = requests.get(link, timeout=5)
        response.raise_for_status()  # Will raise an exception for HTTP errors
    except requests.exceptions.Timeout:
        # print("The request timed out")
        display_error("The request timed out")
        return []
    except requests.exceptions.RequestException as err:
        # print("Something went wrong: %s", err)
        display_error(f"Something went wrong: {err}")
        return []

    try:
        data = response.json()
    except json.JSONDecodeError:
        display_error("Invalid JSON response")
        # print("Invalid JSON response")
        return []

    if "statistics" not in data:
        display_error("Missing 'statistics' in response")
        # print("Missing 'statistics' in response")
        return []

    bases = [entry["libraryName"] for entry in data["statistics"]]
    bases.sort()
    return bases


def run_rra(ranks):
    rra_lib = importr("RobustRankAggreg")
    r_source = ro.r["source"]
    r_source("workflow/notebooks/functions_R.R")
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_ranks = ro.conversion.py2rpy(ranks)
    RobustRankAggregate = ro.r["RobustRankAggregate"]
    res = RobustRankAggregate(r_ranks)
    print(res)


def genemania_link_results(token, tools_available):
    config = "./config/%s.yaml" % token
    content = open_yaml(config)
    if content["screen_type"] == "not_gene":
        return HTML(
            """<p style="color:red;background-color: white;padding: 0.5em;">This module is available only if genes symbol are available.</p>"""
        )

    def update_conditions(update):
        conditions_list = list(tools_available[tool.value].keys())
        conditions.options = conditions_list
        conditions.value = conditions_list[0]

    def on_button_clicked(event):
        treatment, baseline = conditions.value.split("_vs_")
        print("Baseline:", baseline)
        print("Treatment:", treatment)
        print("FDR cut-off:", fdr_cutoff.value)
        print("Score cut-off", score_cutoff.value)
        print("Tool:", tool.value)

        tool_res = tools_available[tool.value]

        if tool.value == "MAGeCK_MLE":
            info = tool_res[conditions.value][conditions.value + ".gene_summary.txt"]
            info = info.loc[info[treatment + "|fdr"] < fdr_cutoff.value]
            if float(score_cutoff.value) < 0:
                info = info.loc[info[treatment + "|beta"] < score_cutoff.value]
            elif float(score_cutoff.value) > 0:
                info = info.loc[info[treatment + "|beta"] > score_cutoff.value]
            genes = info["Gene"]

        if tool.value == "MAGeCK_RRA":
            info = tool_res[conditions.value][conditions.value + ".gene_summary.txt"]
            info = info.loc[info["neg|fdr"] < fdr_cutoff.value]
            genes = info["id"]

        if tool.value == "BAGEL":
            info = tool_res[conditions.value][conditions.value + "_BAGEL_output.pr"]
            info = info.loc[info["BF"] > score_cutoff.value]
            genes = info["Gene"]

        if tool.value == "directional_scoring_method":
            info = tool_res[conditions.value][
                conditions.value + "_all-elements_directional_scoring_method.txt"
            ]
            if score_cutoff.value > 0:
                info = info.loc[info["score"] > score_cutoff.value]
            else:
                info = info.loc[info["score"] < score_cutoff.value]
            genes = info["Gene"]

        if tool.value == "SSREA":
            info = tool_res[conditions.value][
                conditions.value + "_all-elements_SSREA.txt"
            ]
            if float(score_cutoff.value) > 0:
                info = info.loc[info["NES"] > float(score_cutoff.value)]
            elif float(score_cutoff.value) < 0:
                info = info.loc[info["NES"] < float(score_cutoff.value)]
            info = info.loc[info["padj"] < fdr_cutoff.value]
            genes = info["pathway"]

        print("Size geneset:", len(genes))
        link = "http://genemania.org/search/homo-sapiens/" + "/".join(genes)
        print(link)

    TOOLS = [tool for tool in tools_available.keys() if tool != "DESeq2"]

    tool = widgets.Dropdown(options=TOOLS, value=TOOLS[0], description="Tool:")
    tool.observe(update_conditions, "value")

    conditions_list = list(tools_available[tool.value].keys())
    conditions = widgets.Dropdown(
        options=conditions_list, value=conditions_list[0], description="Condition:"
    )

    fdr_cutoff = widgets.FloatSlider(
        min=0.0, max=1.0, step=0.01, value=0.05, description="FDR cut-off:"
    )
    score_cutoff = widgets.IntText(value=0, description="Score cut-off:")

    button = widgets.Button(description="Genemania!")
    button.on_click(on_button_clicked)

    display(tool, conditions, fdr_cutoff, score_cutoff, button)


def ranking(treatment, control, token, tools_available, params, bagel_version):
    def get_occurence_df(data):
        essential_genes = []
        for key in list(data.keys()):
            essential_genes.extend(data[key])
        df = pd.DataFrame(
            np.zeros((len(set(essential_genes)), len(data.keys()))),
            index=set(essential_genes),
            columns=data.keys(),
            dtype=int,
        )
        for tool in data.keys():
            for gene in set(essential_genes):
                if gene in data[tool]:
                    df.loc[gene, tool] = 1
        return df

    config = "./config/%s.yaml" % token
    content = open_yaml(config)

    tool_results = {}
    pdList = []
    tool_genes = []

    comparison = treatment + "_vs_" + control

    if params["MAGeCK_MLE"]["on"]:
        score = params["MAGeCK_MLE"]["score"]
        fdr = params["MAGeCK_MLE"]["fdr"]
        greater = params["MAGeCK_MLE"]["greater"]
        mle = tools_available["MAGeCK_MLE"][comparison][
            comparison + ".gene_summary.txt"
        ]
        if not greater:
            mle = mle[
                (mle["%s|beta" % treatment] < score) & (mle["%s|fdr" % treatment] < fdr)
            ]
        else:
            mle = mle[
                (mle["%s|beta" % treatment] > score) & (mle["%s|fdr" % treatment] < fdr)
            ]
        mle["default_rank"] = mle[treatment + "|beta"].rank(method="dense").copy()
        mle = mle[["Gene", "default_rank"]].rename(
            columns={"Gene": "id", "default_rank": "mle_rank"}
        )
        mle_genes = list(mle.id)
        tool_results["MAGeCK MLE"] = mle_genes
        tool_genes.append(mle_genes)

    if params["MAGeCK_RRA"]["on"]:
        score = params["MAGeCK_RRA"]["score"]
        fdr = params["MAGeCK_RRA"]["fdr"]
        greater = params["MAGeCK_RRA"]["greater"]
        direction = params["MAGeCK_RRA"]["direction"]
        rra = tools_available["MAGeCK_RRA"][comparison][
            comparison + ".gene_summary.txt"
        ]
        if direction == "Negative":
            if not greater:
                rra = rra[(rra["neg|lfc"] < score) & (rra["neg|fdr"] < fdr)]
            else:
                rra = rra[(rra["neg|lfc"] > score) & (rra["neg|fdr"] < fdr)]
            rra = rra[["id", "neg|rank"]].rename(columns={"neg|rank": "rra_rank"})
        else:
            if not greater:
                rra = rra[(rra["pos|lfc"] < score) & (rra["pos|fdr"] < fdr)]
            else:
                rra = rra[(rra["pos|lfc"] > score) & (rra["pos|fdr"] < fdr)]
            rra = rra[["id", "pos|rank"]].rename(columns={"pos|rank": "rra_rank"})
        rra_genes = list(rra.id)
        tool_results["MAGeCK RRA"] = rra_genes
        tool_genes.append(rra_genes)

    if params["BAGEL"]["on"]:
        if bagel_version == 2:
            end_file = "_BAGEL_output.pr"
        elif bagel_version == 1:
            end_file = "_BAGEL1_output.bf"
        score = params["BAGEL"]["score"]
        bagel = tools_available["BAGEL"][comparison][comparison + end_file]

        # Rename GENE column to Gene
        bagel = bagel.rename(columns={"GENE": "Gene"})

        bagel = bagel[(bagel["BF"] > score)]
        bagel["default_rank"] = bagel["BF"].rank(method="dense", ascending=False).copy()
        bagel = bagel[["Gene", "default_rank"]].rename(
            columns={"Gene": "id", "default_rank": "bagel_rank"}
        )
        bagel_genes = list(bagel.id)
        tool_results["BAGEL"] = bagel_genes
        tool_genes.append(bagel_genes)

    if params["directional_scoring_method"]["on"]:
        directional_scoring_method = tools_available["directional_scoring_method"][
            comparison
        ][comparison + "_all-elements_directional_scoring_method.txt"]
        if params["directional_scoring_method"]["direction"] == "Negative":
            directional_scoring_method = directional_scoring_method.loc[
                directional_scoring_method.category == "down"
            ]
        elif params["directional_scoring_method"]["direction"] == "Positive":
            directional_scoring_method = directional_scoring_method.loc[
                directional_scoring_method.category == "up"
            ]
        directional_scoring_method["default_rank"] = (
            directional_scoring_method["score"].rank(method="first").copy()
        )
        directional_scoring_method = directional_scoring_method[
            ["Gene", "default_rank"]
        ].rename(
            columns={"Gene": "id", "default_rank": "directional_scoring_method_rank"}
        )
        directional_scoring_method_genes = list(directional_scoring_method.id)
        tool_results["Directional Scoring Method"] = directional_scoring_method_genes
        tool_genes.append(directional_scoring_method_genes)

    if params["SSREA_like"]["on"]:
        score = params["SSREA_like"]["score"]
        fdr = params["SSREA_like"]["fdr"]
        greater = params["SSREA_like"]["greater"]
        ssrea = tools_available["SSREA"][comparison][
            comparison + "_all-elements_SSREA.txt"
        ]
        if not greater:
            ssrea = ssrea[(ssrea["NES"] < score) & (ssrea["padj"] < fdr)]
        else:
            ssrea = ssrea[(ssrea["NES"] > score) & (ssrea["padj"] < fdr)]
        ssrea["default_rank"] = ssrea["NES"].rank(method="dense").copy()
        ssrea = ssrea[["pathway", "default_rank"]].rename(
            columns={"pathway": "id", "default_rank": "ssrea_rank"}
        )
        ssrea_genes = list(ssrea.id)
        tool_results["SSREA"] = ssrea_genes
        tool_genes.append(ssrea_genes)

    if params["CRISPhieRmix"]["on"]:
        score = params["CRISPhieRmix"]["mean_log2FoldChange"]
        fdr = params["CRISPhieRmix"]["fdr"]
        greater = params["CRISPhieRmix"]["greater"]
        crisphie = tools_available["CRISPhieRmix"][comparison][comparison + ".txt"]
        if not greater:
            crisphie = crisphie[
                (crisphie["mean_log2FoldChange"] < score) & (crisphie["locfdr"] < fdr)
            ]
        else:
            crisphie = crisphie[
                (crisphie["mean_log2FoldChange"] > score) & (crisphie["locfdr"] < fdr)
            ]
        crisphie["default_rank"] = crisphie["locfdr"].rank(method="dense").copy()
        crisphie = crisphie[["gene", "default_rank"]].rename(
            columns={"gene": "id", "default_rank": "crisphiermix_rank"}
        )
        crisphie_genes = list(crisphie.id)
        tool_results["CRISPhieRmix"] = crisphie_genes
        tool_genes.append(crisphie_genes)

    l = []
    for genes in tool_genes:
        for gene in genes:
            l.append(gene)

    if params["MAGeCK_MLE"]["on"]:
        mle = tools_available["MAGeCK_MLE"][comparison][
            comparison + ".gene_summary.txt"
        ]
        mle["default_rank"] = mle[treatment + "|beta"].rank(method="dense").copy()
        mle = mle[["Gene", "default_rank"]].rename(
            columns={"Gene": "id", "default_rank": "mle_rank"}
        )
        pdList.append(mle)

    if params["MAGeCK_RRA"]["on"]:
        rra = tools_available["MAGeCK_RRA"][comparison][
            comparison + ".gene_summary.txt"
        ]
        if params["MAGeCK_RRA"]["direction"] == "Negative":
            rra = rra[["id", "neg|rank"]].rename(columns={"neg|rank": "rra_rank"})
        elif params["MAGeCK_RRA"]["direction"] == "Positive":
            rra = rra[["id", "pos|rank"]].rename(columns={"pos|rank": "rra_rank"})
        pdList.append(rra)

    if params["BAGEL"]["on"]:
        bagel = tools_available["BAGEL"][comparison][comparison + end_file]
        # Renew GENE column to Gene
        bagel = bagel.rename(columns={"GENE": "Gene"})
        bagel["default_rank"] = bagel["BF"].rank(method="dense", ascending=False).copy()
        bagel = bagel[["Gene", "default_rank"]].rename(
            columns={"Gene": "id", "default_rank": "bagel_rank"}
        )
        pdList.append(bagel)

    if params["directional_scoring_method"]["on"]:
        directional_scoring_method = tools_available["directional_scoring_method"][
            comparison
        ][comparison + "_all-elements_directional_scoring_method.txt"]
        directional_scoring_method["default_rank"] = (
            directional_scoring_method["score"].rank(method="first").copy()
        )
        directional_scoring_method = directional_scoring_method[
            ["Gene", "default_rank"]
        ].rename(
            columns={"Gene": "id", "default_rank": "directional_scoring_method_rank"}
        )
        pdList.append(directional_scoring_method)

    if params["SSREA_like"]["on"]:
        ssrea = tools_available["SSREA"][comparison][
            comparison + "_all-elements_SSREA.txt"
        ]
        ssrea["default_rank"] = ssrea["NES"].rank(method="dense").copy()
        ssrea = ssrea[["pathway", "default_rank"]].rename(
            columns={"pathway": "id", "default_rank": "ssrea_rank"}
        )
        pdList.append(ssrea)

    if params["CRISPhieRmix"]["on"]:
        crisphie = tools_available["CRISPhieRmix"][comparison][comparison + ".txt"]
        crisphie["default_rank"] = crisphie["locfdr"].rank(method="dense").copy()
        crisphie = crisphie[["gene", "default_rank"]].rename(
            columns={"gene": "id", "default_rank": "crisphiermix_rank"}
        )
        pdList.append(crisphie)

    df_merged_reduced = reduce(
        lambda left, right: pd.merge(left, right, on=["id"], how="outer"), pdList
    )

    ranks = df_merged_reduced[df_merged_reduced["id"].isin(l)]
    occurences = get_occurence_df(tool_results)

    return (ranks, occurences)


def plot_venn(occurences):
    """Plot a venn diagram from the occurences dataframe"""
    venn_lib = importr("venn")
    grdevices = importr("grDevices")
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_occurences = ro.conversion.py2rpy(occurences)

    with rpy2.robjects.lib.grdevices.render_to_bytesio(
        grdevices.png, width=4096, height=3584, res=600
    ) as img:
        venn_lib.venn(
            r_occurences,
            ilabels=False,
            zcolor="style",
            ilcs=1,
            sncs=1,
            borders=False,
            box=False,
        )
    display(
        IPython.display.display(
            IPython.display.Image(data=img.getvalue(), format="png", embed=True)
        )
    )


def plot_upset(occurences):
    """Plot an upset plot from the occurences dataframe"""
    # Convert 0 and 1 to False and True
    occurences = occurences.astype(bool)

    occurences = occurences.groupby(occurences.columns.tolist()).size()

    # Create the upset plot
    upset = upsetplot.UpSet(occurences)

    # Display the upset plot in the notebook
    upset.plot()
    plt.show()


def reset_params():
    params = {
        "MAGeCK_MLE": {"on": False, "fdr": 0.05, "score": 0, "greater": False},
        "MAGeCK_RRA": {
            "on": False,
            "fdr": 0.05,
            "score": 0,
            "greater": False,
            "direction": "Negative",
        },
        "BAGEL": {"on": False, "score": 0, "greater": True},
        "CRISPhieRmix": {
            "on": False,
            "fdr": 0.05,
            "mean_log2FoldChange": 0,
            "greater": False,
        },
        "directional_scoring_method": {"on": False, "direction": "Negative"},
        "SSREA_like": {"on": False, "fdr": 0.25, "score": 0, "greater": False},
    }
    return params


def display_tools_widgets(tools_selected):
    ### MLE
    def mle_order_update(change):
        if change["new"] == "Greater than score":
            params["MAGeCK_MLE"]["greater"] = True
        else:
            params["MAGeCK_MLE"]["greater"] = False

    def mle_fdr_update(change):
        params["MAGeCK_MLE"]["fdr"] = change["new"]

    def mle_score_update(change):
        params["MAGeCK_MLE"]["score"] = change["new"]

    ### RRA
    def rra_order_update(change):
        if change["new"] == "Greater than score":
            params["MAGeCK_RRA"]["greater"] = True
        else:
            params["MAGeCK_RRA"]["greater"] = False

    def rra_fdr_update(change):
        params["MAGeCK_RRA"]["fdr"] = change["new"]

    def rra_score_update(change):
        params["MAGeCK_RRA"]["score"] = change["new"]

    def rra_direction_update(change):
        params["MAGeCK_RRA"]["direction"] = change["new"]

    ### BAGEL
    def bagel_order_update(change):
        if change["new"] == "Greater than score":
            params["BAGEL"]["greater"] = True
        else:
            params["BAGEL"]["greater"] = False

    def bagel_fdr_update(change):
        params["BAGEL"]["fdr"] = change["new"]

    def bagel_score_update(change):
        params["BAGEL"]["score"] = change["new"]

    ### CRISPhieRmix
    def CRISPhieRmix_order_update(change):
        if change["new"] == "Greater than score":
            params["CRISPhieRmix"]["greater"] = True
        else:
            params["CRISPhieRmix"]["greater"] = False

    def CRISPhieRmix_fdr_update(change):
        params["CRISPhieRmix"]["fdr"] = change["new"]

    def CRISPhieRmix_score_update(change):
        params["CRISPhieRmix"]["mean_log2FoldChange"] = change["new"]

    ### directional_scoring_method
    def directional_scoring_method_direction_update(change):
        params["directional_scoring_method"]["direction"] = change["new"]

    ### SSREA
    def SSREA_like_order_update(change):
        if change["new"] == "Greater than score":
            params["SSREA_like"]["greater"] = True
        else:
            params["SSREA_like"]["greater"] = False

    def SSREA_like_fdr_update(change):
        params["SSREA_like"]["fdr"] = change["new"]

    def SSREA_like_score_update(change):
        params["SSREA_like"]["score"] = change["new"]

    if "MAGeCK_MLE" in tools_selected:
        params["MAGeCK_MLE"]["on"] = True
        mle_fdr = widgets.FloatSlider(
            min=0.0, max=1.0, step=0.01, value=0.05, description="FDR:"
        )
        mle_score = widgets.FloatText(value=0, description="Score cut-off:")
        mle_text = widgets.HTML(value="<b>MAGeCK MLE</b>:")
        mle_order = widgets.ToggleButtons(
            options=["Lower than score", "Greater than score"],
            description="Selection:",
            name="test",
        )
        display(mle_text)
        mle_box = widgets.HBox([mle_fdr, mle_score, mle_order])
        display(mle_box)
        mle_order.observe(mle_order_update, "value")
        mle_score.observe(mle_score_update, "value")
        mle_fdr.observe(mle_fdr_update, "value")
    if "MAGeCK_RRA" in tools_selected:
        params["MAGeCK_RRA"]["on"] = True
        rra_direction = widgets.ToggleButtons(
            options=["Negative", "Positive"], description="Direction:"
        )
        rra_fdr = widgets.FloatSlider(
            min=0.0, max=1.0, step=0.01, value=0.05, description="FDR:"
        )
        rra_score = widgets.FloatText(value=0, description="Score cut-off:")
        rra_text = widgets.HTML(value="<b>MAGeCK RRA</b>:")
        rra_order = widgets.ToggleButtons(
            options=["Lower than score", "Greater than score"], description="Selection:"
        )
        display(rra_text)
        rra_box = widgets.HBox([rra_fdr, rra_direction, rra_score, rra_order])
        display(rra_box)
        rra_direction.observe(rra_direction_update, "value")
        rra_order.observe(rra_order_update, "value")
        rra_score.observe(rra_score_update, "value")
        rra_fdr.observe(rra_fdr_update, "value")
    if "BAGEL" in tools_selected:
        params["BAGEL"]["on"] = True
        bagel_score = widgets.FloatText(
            value=0,
            description="BF >",
            layout=layout,
            display="flex",
            flex_flow="column",
            align_items="stretch",
        )
        bagel_text = widgets.HTML(value="<b>BAGEL</b>:")
        # bagel_order = widgets.ToggleButtons(
        #     options=["Greater than score", "Lower than score"], description="Selection:"
        # )
        display(bagel_text)
        bagel_box = widgets.HBox([bagel_score])  # , bagel_order])
        display(bagel_box)
        # bagel_order.observe(bagel_order_update, "value")
        bagel_score.observe(bagel_score_update, "value")
    if "CRISPhieRmix" in tools_selected:
        params["CRISPhieRmix"]["on"] = True
        CRISPhieRmix_fdr = widgets.FloatSlider(
            min=0.0, max=1.0, step=0.01, value=0.05, description="FDR:"
        )
        CRISPhieRmix_score = widgets.FloatText(value=0, description="Score cut-off:")
        CRISPhieRmix_text = widgets.HTML(value="<b>CRISPhieRmix</b>:")
        CRISPhieRmix_order = widgets.ToggleButtons(
            options=["Lower than score", "Greater than score"], description="Selection:"
        )
        display(CRISPhieRmix_text)
        CRISPhieRmix_box = widgets.HBox(
            [CRISPhieRmix_fdr, CRISPhieRmix_score, CRISPhieRmix_order]
        )
        display(CRISPhieRmix_box)
        CRISPhieRmix_order.observe(CRISPhieRmix_order_update, "value")
        CRISPhieRmix_score.observe(CRISPhieRmix_score_update, "value")
        CRISPhieRmix_fdr.observe(CRISPhieRmix_fdr_update, "value")
    if "directional_scoring_method" in tools_selected:
        params["directional_scoring_method"]["on"] = True
        directional_scoring_method_direction = widgets.ToggleButtons(
            options=["Negative", "Positive"], description="Direction:"
        )
        directional_scoring_method_text = widgets.HTML(
            value="<b>Directional Scoring Method</b>:"
        )
        display(directional_scoring_method_text)
        directional_scoring_method_box = widgets.HBox(
            [directional_scoring_method_direction]
        )
        display(directional_scoring_method_box)
        directional_scoring_method_direction.observe(
            directional_scoring_method_direction_update, "value"
        )
    if "SSREA" in tools_selected:
        params["SSREA_like"]["on"] = True
        SSREA_like_fdr = widgets.FloatSlider(
            min=0.0, max=1.0, step=0.01, value=0.25, description="FDR:"
        )
        SSREA_like_score = widgets.FloatText(value=0, description="Score cut-off:")
        SSREA_like_text = widgets.HTML(value="<b>SSREA</b>:")
        SSREA_like_order = widgets.ToggleButtons(
            options=["Lower than score", "Greater than score"], description="Selection:"
        )
        display(SSREA_like_text)
        SSREA_like_box = widgets.HBox(
            [SSREA_like_fdr, SSREA_like_score, SSREA_like_order]
        )
        display(SSREA_like_box)
        SSREA_like_order.observe(SSREA_like_order_update, "value")
        SSREA_like_score.observe(SSREA_like_score_update, "value")
        SSREA_like_fdr.observe(SSREA_like_fdr_update, "value")


def disable_widgets(token):
    config = "./config/%s.yaml" % token
    content = open_yaml(config)
    disabled = False
    if content["screen_type"] == "not_gene" and content["bed_annotation_file"] == "":
        disabled = True
    return disabled


def multiple_tools_results(tools_available, token, bagel_version):
    TOOLS = [tool for tool in tools_available.keys() if not tool in ["DESeq2"]]

    # Define widgets's options
    conditions_options = tools_available[TOOLS[0]].keys()
    tools_options = TOOLS

    # Define widgets
    conditions_widget = widgets.Dropdown(
        options=conditions_options, description="Conditions:"
    )
    tools_widget = widgets.SelectMultiple(
        options=tools_options,
        description="Tool:",
    )
    selection_widgets = widgets.ToggleButtons(
        options=["Intersection", "Union"],
        description="Selection mode:",
        tooltips=[
            "Use elements at intersection of all selected methods",
            "Use union of elements of all selected methods",
        ],
    )

    output_tools_form = widgets.Output()

    def plot_interactive_heatmap(table, column_to_plot):
        """Plot an interactive heatmap with dendrogram using Plotly from a depmap table."""

        # Transform datatable to a wide form.
        data = (
            pd.pivot_table(
                table,
                index="cell_line_name",
                columns="gene_name",
                values=column_to_plot,
            )
            .reset_index()
            .set_index("cell_line_name")
            .dropna(axis=1)
        )

        # Store datable index and columns in two lists
        labels = data.index
        cols = data.columns

        # Convert datatable to numpy array
        data = data.to_numpy(dtype=float)

        # Compute z-score.
        if column_to_plot in ["rna_expression", "protein_expression"]:
            data_array = stats.zscore(data, axis=0)
        else:
            data_array = data

        # Initialize figure by creating upper dendrogram
        fig = ff.create_dendrogram(data_array, orientation="bottom")
        for i in range(len(fig["data"])):
            fig["data"][i]["yaxis"] = "y2"

        # Create Side Dendrogram
        dendro_side = ff.create_dendrogram(
            data_array.transpose(),
            orientation="right",
        )  # , labels=cols)
        for i in range(len(dendro_side["data"])):
            dendro_side["data"][i]["xaxis"] = "x2"

        # Add Side Dendrogram Data to Figure
        for data in dendro_side["data"]:
            fig.add_trace(data)

        # Create Heatmap
        dendro_leaves_side = dendro_side["layout"]["yaxis"]["ticktext"]
        dendro_leaves_side = list(map(int, dendro_leaves_side))

        dendro_leaves_upper = fig["layout"]["xaxis"]["ticktext"]
        dendro_leaves_upper = list(map(int, dendro_leaves_upper))

        heat_data = data_array.transpose()
        heat_data = heat_data[dendro_leaves_side, :]
        heat_data = heat_data[:, dendro_leaves_upper]

        # Center colorbar on zero
        max_value = np.max(heat_data)
        min_value = np.min(heat_data)
        prop_zero = abs(min_value / (abs(max_value) + abs(min_value)))

        # Define heatmap
        heatmap = [
            go.Heatmap(
                x=dendro_leaves_upper,
                y=dendro_leaves_side,
                z=heat_data,
                colorscale=[
                    [0.0, "rgb(0, 13, 255)"],
                    [prop_zero, "rgb(255, 255, 255)"],
                    [1.0, "rgb(255, 0, 38)"],
                ],
                colorbar={
                    "x": -0.1,
                    "tickfont": {"size": 14},
                    "lenmode": "fraction",
                    "len": 0.25,
                    "thickness": 20,
                },
            )
        ]

        heatmap[0]["x"] = fig["layout"]["xaxis"]["tickvals"]
        heatmap[0]["y"] = dendro_side["layout"]["yaxis"]["tickvals"]

        indicies_xaxis = list(map(int, fig["layout"]["xaxis"]["ticktext"]))
        indicies_yaxis = list(map(int, dendro_side["layout"]["yaxis"]["ticktext"]))

        # Add Heatmap Data to Figure
        for data in heatmap:
            fig.add_trace(data)

        # Edit Layout
        fig.update_layout(
            {"width": 800, "height": 900, "showlegend": False, "hovermode": "closest"}
        )

        fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)")

        # Edit xaxis
        fig.update_layout(
            font={"size": 10},
            xaxis={
                "domain": [0.15, 1],
                "ticktext": labels[indicies_xaxis],
                "mirror": False,
                "showgrid": False,
                "showline": False,
                "zeroline": False,
                "showticklabels": True,
                "tickangle": 45,
                "ticks": "",
            },
        )
        # Edit xaxis2
        fig.update_layout(
            xaxis2={
                "domain": [0, 0.15],
                "mirror": False,
                "showgrid": False,
                "showline": False,
                "zeroline": False,
                "showticklabels": False,
                "ticks": "",
            }
        )

        # Wrap elements in cols[indicies_yaxis] with <i>
        wrapped_cols = pd.Index([f"<i>{col}</i>" for col in cols], name=cols.name)

        # # Edit yaxis
        fig.update_layout(
            font={"size": 10},
            yaxis={
                "domain": [0, 0.85],
                "ticktext": wrapped_cols[indicies_yaxis],
                "tickmode": "array",
                "tickvals": dendro_side["layout"]["yaxis"]["tickvals"],
                "mirror": False,
                "showgrid": False,
                "showline": False,
                "zeroline": False,
                "showticklabels": True,
                "ticks": "",
                "side": "right",
            },
        )
        # Edit yaxis2
        fig.update_layout(
            yaxis2={
                "domain": [0.825, 0.975],
                "mirror": False,
                "showgrid": False,
                "showline": False,
                "zeroline": False,
                "showticklabels": False,
                "ticks": "",
            }
        )

        # Plot!
        display(
            HTML(
                """<h4 style="color:white;font-weight: bold;background-color: red;padding: 0.5em;">Depmap %s heatmap</h4>"""
                % column_to_plot
            )
        )
        show_parameters(params)
        display(fig)

        # Create download button with plotly.io.write_image
        download_button = widgets.Button(
            description="Download", style=dict(button_color="#707070")
        )

        @download_button.on_click
        def download_button_clicked(b):
            heatmap_name = f"results/{token}/{column_to_plot}_heatmap.pdf"
            fig.write_image(heatmap_name)
            # Get the base name of the heatmap
            base_name = os.path.basename(heatmap_name)

            # Display a clickable link to the image
            display(
                HTML(
                    f'<a href="{base_name}" target="_blank">Open {column_to_plot} heatmap</a>'
                )
            )

        display(download_button)

    @output_tools_form.capture()
    def tools_widget_updated(event):
        global params
        params = reset_params()

        output_tools_form.clear_output()
        display_tools_widgets(event["new"])

    tools_widget.observe(tools_widget_updated, "value")
    display(widgets.VBox([conditions_widget, tools_widget, selection_widgets]))
    display(output_tools_form)

    disabled = disable_widgets(token)
    venn_button = widgets.Button(
        description="Venn/Upset plots", style=dict(button_color="#707070")
    )
    rra_button = widgets.Button(
        description="RRA ranking", style=dict(button_color="#707070")
    )
    genemania_button = widgets.Button(
        description="Genemania", style=dict(button_color="#707070"), disabled=disabled
    )
    enrichr_button = widgets.Button(
        description="EnrichR", style=dict(button_color="#707070"), disabled=disabled
    )
    depmap_button = widgets.Button(
        description="depmap", style=dict(button_color="#707070"), disabled=disabled
    )
    ranking_button = widgets.Button(
        description="Save  ranking", style=dict(button_color="#707070")
    )

    buttons_box = widgets.HBox(
        [
            venn_button,
            rra_button,
            genemania_button,
            enrichr_button,
            depmap_button,
            ranking_button,
        ]
    )
    display(buttons_box)

    def show_parameters(params):
        if "MAGeCK_MLE" in tools_widget.value:
            if params["MAGeCK_MLE"]["greater"]:
                word = "greater"
            else:
                word = "less"
            display(
                HTML(
                    """<p style="color:black;padding: 0.5em;"><b>MAGeCK MLE</b>: FDR < %s, β score threshold = %s, keep elements with β score %s than threshold.</p>"""
                    % (params["MAGeCK_MLE"]["fdr"], params["MAGeCK_MLE"]["score"], word)
                )
            )
        if "MAGeCK_RRA" in tools_widget.value:
            if params["MAGeCK_RRA"]["greater"]:
                word = "greater"
            else:
                word = "less"
            display(
                HTML(
                    """<p style="color:black;padding: 0.5em;"><b>MAGeCK RRA</b>: FDR < %s, LogFoldChange threshold = %s, keep elements with LogFoldChange %s than threshold and use FDR values associated to %s selection.</p>"""
                    % (
                        params["MAGeCK_RRA"]["fdr"],
                        params["MAGeCK_RRA"]["score"],
                        word,
                        params["MAGeCK_RRA"]["direction"],
                    )
                )
            )
        if "BAGEL" in tools_widget.value:
            # if params["BAGEL"]["greater"]:
            #     word = "greater"
            # else:
            #     word = "less"
            display(
                HTML(
                    """<p style="color:black;padding: 0.5em;"><b>BAGEL</b>: Bayesian factor threshold = %s, keep elements with Bayesian factor > than threshold.</p>"""
                    % (params["BAGEL"]["score"])
                )
            )
        if "CRISPhieRmix" in tools_widget.value:
            if params["CRISPhieRmix"]["greater"]:
                word = "greater"
            else:
                word = "less"
            display(
                HTML(
                    """<p style="color:black;padding: 0.5em;"><b>CRISPhieRmix</b>: FDR < %s, Log2FoldChange Average threshold = %s, keep elements with Log2FoldChange Average %s than threshold.</p>"""
                    % (
                        params["CRISPhieRmix"]["fdr"],
                        params["CRISPhieRmix"]["mean_log2FoldChange"],
                        word,
                    )
                )
            )
        if "SSREA" in tools_widget.value:
            if params["SSREA_like"]["greater"]:
                word = "greater"
            else:
                word = "less"
            display(
                HTML(
                    """<p style="color:black;padding: 0.5em;"><b>SSREA (fGSEA)</b>: FDR < %s, NES threshold = %s, keep elements with NES %s than threshold.</p>"""
                    % (params["SSREA_like"]["fdr"], params["SSREA_like"]["score"], word)
                )
            )
        if "directional_scoring_method" in tools_widget.value:
            display(
                HTML(
                    """<p style="color:black;padding: 0.5em;"><b>Directional Scoring Method</b>: %s selection.</p>"""
                    % (params["directional_scoring_method"]["direction"])
                )
            )

    def venn_button_clicked(b):
        # If at least one tool is selected
        if len(tools_widget.value) > 0:
            treatment, control = conditions_widget.value.split("_vs_")
            ranks, occurences = ranking(
                treatment, control, token, tools_available, params, bagel_version
            )
            # if selection_widgets.value == "Intersection":
            df = pd.DataFrame(
                occurences.eq(occurences.iloc[:, 0], axis=0).all(1),
                columns=["intersection"],
            )
            genes_list_at_intersection = df.loc[df.intersection == True].index
            # else:
            df = pd.DataFrame(
                occurences.eq(occurences.iloc[:, 0], axis=0).any(1), columns=["union"]
            )
            genes_list_at_union = df.loc[df.union == True].index
            display(
                HTML(
                    """<p style="color:white;font-weight: bold;background-color: orange;padding: 0.5em;">Venn diagram: %s</p>"""
                    % selection_widgets.value
                )
            )

            show_parameters(params)
            display(occurences.sum(axis=0))
            plot_venn(occurences)
            # If more than one tool is selected, display the upset plot
            if len(tools_widget.value) > 1:
                plot_upset(occurences)

            textarea_intersection = widgets.Textarea(
                value="\n".join(genes_list_at_intersection),
                description=f"List of genes at intersection of all methods (n = {len(genes_list_at_intersection)}):",
                disabled=True,
                # Increase the height of the textarea (default is 6 rows)
                layout=widgets.Layout(height="200px", width="60%"),
                style={"description_width": "400px"},
            )
            display(textarea_intersection)

            textarea_union = widgets.Textarea(
                value="\n".join(genes_list_at_union),
                description=f"List of genes at union of all methods (n = {len(genes_list_at_union)}):",
                disabled=True,
                # Increase the height of the textarea (default is 6 rows)
                layout=widgets.Layout(height="200px", width="60%"),
                style={"description_width": "400px"},
            )

            display(textarea_union)

            # Utilisez du HTML et du JavaScript pour changer la couleur du texte
            display(
                HTML(
                    """
                <style>
                    .widget-textarea textarea:disabled {
                        color: black !important;
                        opacity: 1 !important;
                    }
                </style>
                <script>
                    require(["base/js/namespace"], function(Jupyter) {
                        Jupyter.notebook.events.one("kernel_ready.Kernel", function() {
                            Jupyter.notebook.execute_cells([Jupyter.notebook.get_cell_index(Jupyter.notebook.get_selected_cell())]);
                        });
                    });
                </script>
            """
                )
            )
        else:
            display_warning("Please select at least one tool above.")

    def rra_button_clicked(b):
        if len(tools_widget.value) == 0:
            display_warning("Please select at least one tool above.")
        else:

            treatment, control = conditions_widget.value.split("_vs_")
            ranks, occurences = ranking(
                treatment, control, token, tools_available, params, bagel_version
            )
            display(
                HTML(
                    """<p style="color:white;font-weight: bold;background-color: green;padding: 0.5em;">RRA results</p>"""
                )
            )
            show_parameters(params)
            run_rra(ranks)

    def genemania_button_clicked(b):
        if len(tools_widget.value) == 0:
            display_warning("Please select at least one tool above.")
        else:
            treatment, control = conditions_widget.value.split("_vs_")
            ranks, occurences = ranking(
                treatment, control, token, tools_available, params, bagel_version
            )
            if selection_widgets.value == "Intersection":
                df = pd.DataFrame(
                    occurences.eq(occurences.iloc[:, 0], axis=0).all(1),
                    columns=["intersection"],
                )
                genes_list = df.loc[df.intersection == True].index
            else:
                df = pd.DataFrame(
                    occurences.eq(occurences.iloc[:, 0], axis=0).any(1),
                    columns=["union"],
                )
                genes_list = df.loc[df.union == True].index
            genes_list = regions2genes(token, genes_list)
            display(
                HTML(
                    """<p style="color:white;font-weight: bold;background-color: blue;padding: 0.5em;">Genemania link - %s</p>"""
                    % selection_widgets.value
                )
            )
            show_parameters(params)
            link = "http://genemania.org/search/homo-sapiens/" + "/".join(genes_list)
            print("Link to Genemania website: (%s elements)\n" % len(genes_list))
            print(link)

    def enrichr_button_clicked(b):
        def show_enrichr_plots(
            b,
            genes,
            bases,
            size,
            plot_type,
            col_2,
            col_1,
            description,
            # output
        ):
            charts = []
            title = description.value + " (%s)" % selection_widgets.value
            # with output:
            show_parameters(params)
            for base in bases.value:
                enrichr_res = get_enrichr_results(genes, description.value, base)
                table = create_enrichr_table(enrichr_res)
                if plot_type.value == "Bar":
                    chart = enrichmentBarPlot(
                        table, size.value, title, col_1.value, col_2.value, base
                    )
                else:
                    chart = enrichmentCirclePlot(
                        table, size.value, title, col_1.value, col_2.value, base
                    )
                charts.append(chart)
            for chart in charts:
                display(chart)

        if len(tools_widget.value) == 0:
            display_warning("Please select at least one tool above.")
        else:
            display(
                HTML(
                    f"""
                    <p style="color:white;
                            font-weight: bold;
                            background-color: purple;
                            padding: 0.5em;">EnrichR for genes at {selection_widgets.value}</p>
                    """
                )
            )
            treatment, control = conditions_widget.value.split("_vs_")
            ranks, occurences = ranking(
                treatment, control, token, tools_available, params, bagel_version
            )
            if selection_widgets.value == "Intersection":
                df = pd.DataFrame(
                    occurences.eq(occurences.iloc[:, 0], axis=0).all(1),
                    columns=["intersection"],
                )
                genes_list = df.loc[df.intersection == True].index
            else:
                df = pd.DataFrame(
                    occurences.eq(occurences.iloc[:, 0], axis=0).any(1),
                    columns=["union"],
                )
                genes_list = df.loc[df.union == True].index
            genes_list = regions2genes(token, genes_list)
            BASES = get_enrichr_bases()
            bases = widgets.SelectMultiple(
                options=BASES, description="Genesets:", rows=10
            )
            col_2 = widgets.ColorPicker(
                concise=False, description="Top color", value="blue"
            )
            col_1 = widgets.ColorPicker(
                concise=False, description="Bottom color", value="red"
            )
            plot_type = widgets.Dropdown(
                options=["Circle", "Bar"], value="Circle", description="Plot type:"
            )
            size = widgets.Dropdown(
                options=[5, 10, 20, 50, 100, 200, "max"], value=5, description="Size:"
            )
            description = widgets.Text(
                value="My gene list",
                placeholder="Description",
                description="Description:",
            )
            button_enrichr = widgets.Button(description="Plot!")

            display(
                widgets.VBox(
                    [description, bases, plot_type, size, col_2, col_1, button_enrichr]
                )
            )

            # output_enrichr = widgets.Output()

            # display(output_enrichr)

            button_enrichr.on_click(
                partial(
                    show_enrichr_plots,
                    genes=genes_list,
                    bases=bases,
                    size=size,
                    plot_type=plot_type,
                    col_2=col_2,
                    col_1=col_1,
                    description=description,
                    # output=output_enrichr,
                )
            )

    def depmap_button_clicked(b):
        if len(tools_widget.value) == 0:
            display_warning("Please select at least one tool above.")
        else:
            display(
                HTML(
                    """<p style="color:white;font-weight: bold;background-color: #A52A2A;padding: 0.5em;">DepMap vizualisation module</p>"""
                )
            )
            data_types_widget = widgets.RadioButtons(
                options=["crispr", "proteomic", "rnai", "tpm", "mutations"],
                value="crispr",
                description="Choose a data type:",
                disabled=False,
            )

            def download_depmap_file(data_type, release):
                """Determine if the file is already downloaded or not."""
                target_file = f"resources/depmap/{release}_{data_type}.txt"
                for file_name in listdir("resources/depmap/"):
                    if ("metadata" in file_name) and (not release in file_name):
                        os.remove(file_name)
                    if file_name == target_file:
                        return False
                    elif data_type in file_name:
                        if not release in str(file_name):
                            os.remove(file_name)
                            return True
                        else:
                            return False
                return True

            def get_release():
                """Get the release of the depmap data."""
                depmap_release = depmap.depmap_release()
                return str(depmap_release).rstrip()[5:-1]

            def depmap_query_button_clicked(b):
                """Query the depmap data."""
                # Get the release of the depmap data
                depmap_release = get_release()
                # Set the path to save the file
                save_path = (
                    f"resources/depmap/{depmap_release}_{data_types_widget.value}.txt"
                )
                # Create the directory if it does not exist
                Path("resources/depmap/").mkdir(parents=True, exist_ok=True)
                # Get the treatment and control conditions
                treatment, control = conditions_widget.value.split("_vs_")
                # Get the ranking
                ranks, occurences = ranking(
                    treatment, control, token, tools_available, params, bagel_version
                )
                if download_depmap_file(data_types_widget.value, depmap_release):
                    # print("This step can take some time.")
                    display_info("This step can take some time.")
                    # print("Querying: %s..." % data_types_widget.value)
                    display_info(f"Querying: {data_types_widget.value}...", bold=False)
                    eh = experimentHub.ExperimentHub()
                    base_package = importr("base")
                    dplyr = importr("dplyr")
                    tidyr = importr("tidyr")
                    readr = importr("readr")
                    if data_types_widget.value == "rnai":
                        depmap_data = tidyr.drop_na(
                            dplyr.select(
                                depmap.depmap_rnai(),
                                "depmap_id",
                                "gene_name",
                                "dependency",
                            )
                        )
                    elif data_types_widget.value == "crispr":
                        depmap_data = tidyr.drop_na(
                            dplyr.select(
                                depmap.depmap_crispr(),
                                "depmap_id",
                                "gene_name",
                                "dependency",
                            )
                        )
                    elif data_types_widget.value == "proteomic":
                        depmap_data = tidyr.drop_na(
                            dplyr.select(
                                depmap.depmap_proteomic(),
                                "depmap_id",
                                "gene_name",
                                "protein_expression",
                            )
                        )
                    elif data_types_widget.value == "tpm":
                        depmap_data = tidyr.drop_na(
                            dplyr.select(
                                depmap.depmap_TPM(),
                                "depmap_id",
                                "gene_name",
                                "rna_expression",
                            )
                        )
                    elif data_types_widget.value == "mutations":
                        depmap_data = tidyr.drop_na(
                            dplyr.select(
                                depmap.depmap_mutationCalls(),
                                "depmap_id",
                                "gene_name",
                                "protein_change",
                                "is_deleterious",
                            )
                        )
                    if not os.path.isfile(
                        f"resources/depmap/{depmap_release}_metadata.txt"
                    ):
                        depmap_metadata = dplyr.select(
                            depmap.depmap_metadata(),
                            "depmap_id",
                            "sample_collection_site",
                            "primary_or_metastasis",
                            "primary_disease",
                            "subtype_disease",
                            "cell_line",
                            "cell_line_name",
                        )
                        print("Saving metadata...")
                        utils.write_table(
                            depmap_metadata,
                            "resources/depmap/{depmap_release}_metadata.txt",
                            row_names=False,
                            quote=False,
                            sep="\t",
                        )
                    else:
                        print("Import metadata...")
                        depmap_metadata = readr.read_delim(
                            f"resources/depmap/{depmap_release}_metadata.txt",
                            delim="\t",
                        )
                    depmap_data = base_package.merge(
                        depmap_data, depmap_metadata, by="depmap_id"
                    )
                    print(f"Saving {save_path}")
                    utils.write_table(
                        depmap_data, save_path, row_names=False, quote=False, sep="\t"
                    )
                print(f"Opening {save_path}")
                data = pd.read_table(save_path, sep="\t")

                tissues_init = list(set(data.cell_line))
                tissues = [
                    "_".join(str(tissu).split("_")[1:])
                    for tissu in tissues_init
                    if str(tissu) not in ["nan", ""]
                ]
                tissues = list(set(tissues))
                tissues.insert(0, "All")

                primary_diseases = list(set(data.primary_disease))
                primary_diseases.insert(0, "All")

                cell_lines_init = list(set(data.cell_line_name))
                cell_lines = [
                    tissu for tissu in cell_lines_init if str(tissu) not in ["nan", ""]
                ]
                cell_lines = natural_sort(cell_lines)
                cell_lines.insert(0, "All")

                tissues_widget = widgets.SelectMultiple(
                    options=tissues, value=["All"], description="Tissu:"
                )
                primary_diseases_widget = widgets.SelectMultiple(
                    options=primary_diseases,
                    value=["All"],
                    description="Primary tissu:",
                )
                cell_lines_widget = widgets.SelectMultiple(
                    options=cell_lines, value=["All"], description="Cell line:"
                )

                def update_primary_diseases_widget(update):
                    if not "All" in tissues_widget.value:
                        subseted_data = data[
                            data["cell_line"].str.contains(
                                "|".join(tissues_widget.value), na=False
                            )
                        ]
                    else:
                        subseted_data = data
                    primary_diseases = list(set(subseted_data.primary_disease))
                    primary_diseases = natural_sort(primary_diseases)
                    primary_diseases.insert(0, "All")
                    primary_diseases_widget.options = primary_diseases
                    primary_diseases_widget.value = ["All"]

                def update_cell_lines_widget(update):
                    if not "All" in primary_diseases_widget.value:
                        subseted_data = data[
                            data["primary_disease"].str.contains(
                                "|".join(primary_diseases_widget.value), na=False
                            )
                        ]
                    else:
                        if not "All" in tissues_widget.value:
                            subseted_data = data[
                                data["cell_line"].str.contains(
                                    "|".join(tissues_widget.value), na=False
                                )
                            ]
                        else:
                            subseted_data = data
                    cell_lines_init = list(set(subseted_data.cell_line_name))
                    cell_lines = [
                        cell_line
                        for cell_line in cell_lines_init
                        if not str(cell_line) in ["nan"]
                    ]
                    cell_lines = natural_sort(cell_lines)
                    cell_lines.insert(0, "All")
                    cell_lines_widget.options = cell_lines
                    cell_lines_widget.value = ["All"]

                tissues_widget.observe(update_primary_diseases_widget, "value")
                primary_diseases_widget.observe(update_cell_lines_widget, "value")

                def tissu_selection_button_clicked(b):
                    print("Please wait!")
                    dic = {
                        "rnai": "dependency",
                        "crispr": "dependency",
                        "tpm": "rna_expression",
                        "proteomic": "protein_expression",
                        "mutations": ["protein_change", "is_deleterious"],
                    }
                    variable = dic[data_types_widget.value]
                    save_path = "resources/depmap/%s_%s.txt" % (
                        depmap_release,
                        data_types_widget.value,
                    )
                    table = pd.read_table(save_path, sep="\t")
                    if "All" in cell_lines_widget.value:
                        cell_lines_selected = cell_lines_widget.options
                    else:
                        cell_lines_selected = cell_lines_widget.value
                    boolean_series = table.cell_line_name.isin(cell_lines_selected)
                    table = table[boolean_series]
                    if selection_widgets.value == "Intersection":
                        df = pd.DataFrame(
                            occurences.eq(occurences.iloc[:, 0], axis=0).all(1),
                            columns=["intersection"],
                        )
                        genes_list = df.loc[df.intersection == True].index
                    else:
                        df = pd.DataFrame(
                            occurences.eq(occurences.iloc[:, 0], axis=0).any(1),
                            columns=["union"],
                        )
                        genes_list = df.loc[df.union == True].index
                    genes_list = regions2genes(token, genes_list)
                    if data_types_widget.value == "mutations":
                        columns = [
                            "gene_name",
                            "cell_line",
                            "cell_line_name",
                            "sample_collection_site",
                            "primary_or_metastasis",
                            "primary_disease",
                            "subtype_disease",
                        ]
                        for value in variable:
                            columns.append(value)
                        essential_genes = table[table.gene_name.isin(genes_list)][
                            columns
                        ]
                        essential_genes = essential_genes.loc[
                            essential_genes.is_deleterious == True
                        ]
                        essential_genes = (
                            essential_genes.groupby(
                                [
                                    "gene_name",
                                    "cell_line",
                                    "cell_line_name",
                                    "sample_collection_site",
                                    "primary_or_metastasis",
                                    "primary_disease",
                                    "subtype_disease",
                                ]
                            )["protein_change"]
                            .apply(lambda x: "%s" % "|".join(map(str, x)))
                            .reset_index()
                        )
                        chart = (
                            alt.Chart(
                                essential_genes, title="DepMap Deleterious Mutations"
                            )
                            .mark_rect()
                            .encode(
                                x=alt.X(
                                    "cell_line_name", axis=alt.Axis(title="Cell line")
                                ),
                                y=alt.Y("gene_name", axis=alt.Axis(title="Gene name")),
                                color=alt.Color(
                                    "primary_disease",
                                    scale=alt.Scale(scheme="tableau20"),
                                    legend=alt.Legend(title="Primary disease"),
                                ),
                                tooltip=[
                                    "protein_change",
                                    "gene_name",
                                    "cell_line",
                                    "cell_line_name",
                                    "sample_collection_site",
                                    "primary_or_metastasis",
                                    "primary_disease",
                                    "subtype_disease",
                                ],
                            )
                            .interactive()
                        )
                        show_parameters(params)
                        display(chart)
                    else:
                        essential_genes = table[table.gene_name.isin(genes_list)][
                            ["gene_name", "cell_line_name", variable]
                        ]
                        plot_interactive_heatmap(essential_genes, variable)

                button = widgets.Button(description="Run!")
                button.on_click(tissu_selection_button_clicked)
                display(
                    tissues_widget, primary_diseases_widget, cell_lines_widget, button
                )

            depmap_query_button = widgets.Button(description="Query!")
            depmap_query_button.on_click(depmap_query_button_clicked)
            display(data_types_widget, depmap_query_button)

    venn_button.on_click(venn_button_clicked)
    rra_button.on_click(rra_button_clicked)
    genemania_button.on_click(genemania_button_clicked)
    enrichr_button.on_click(enrichr_button_clicked)
    depmap_button.on_click(depmap_button_clicked)

    def ranking_button_clicked(event):
        if len(tools_widget.value) == 0:
            display_warning("Please select at least one tool above.")
        else:
            display(
                HTML(
                    """<p style="color:white;font-weight: bold;background-color: #03fc3d;padding: 0.5em;">Save ranking and occurences</p>"""
                )
            )
            show_parameters(params)
            treatment, control = conditions_widget.value.split("_vs_")
            ranks, occurences = ranking(
                treatment, control, token, tools_available, params, bagel_version
            )
            download_name = widgets.Text(value="analysis", description="Name:")
            download_file(
                content=ranks.to_string(),
                filename="ranking.txt",
                label="Download ranking!",
            )
            download_file(
                content=occurences.to_string(),
                filename="occurences.txt",
                label="Download occurences!",
            )

    ranking_button.on_click(ranking_button_clicked)


def condition_comparison(results_directory, tools_available, token):
    """This function let user choose a tool and two conditions and plot them against each other."""

    def get_data(tool, comparison_1, comparison_2, results_directory):
        """Get the data for the plot."""
        functions = {
            "MAGeCK_MLE": MAGeCK_MLE_data,
            "MAGeCK_RRA": MAGeCK_RRA_data,
            "CRISPhieRmix": CRISPhieRmix_data,
            "SSREA": SSREA_like_data,
            "directional_scoring_method": directional_scoring_method_data,
        }
        data_function = functions[tool]
        comparison_1_data = data_function(
            comparison=comparison_1,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )
        comparison_2_data = data_function(
            comparison=comparison_2,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )
        return comparison_1_data, comparison_2_data

    def parameters_widgets():
        """Create widgets for the parameters."""
        # Define widgets's options
        tools_options = list(tools_available.keys())
        # Remove DESeq2
        if "DESeq2" in tools_options:
            tools_options.remove("DESeq2")
        if "directional_scoring_method" in tools_options:
            tools_options.remove("directional_scoring_method")
        if "SSREA" in tools_options:
            tools_options.remove("SSREA")
        comparisons_options = list(tools_available[tools_options[0]].keys())
        comparisons_options.sort()

        # Create a gene/element selection widget with Text input
        gene_selection_widget = widgets.Text(
            value="",
            placeholder="Comma-separated list",
            description="Gene/Element:",
            disabled=False,
            style={"description_width": "100px"},
        )
        # Create a color widget for the gene/element selection
        color_gene_widget = widgets.ColorPicker(
            concise=True, description="", value="green"
        )

        # Define widgets
        # Create a condition 1 widget
        comparison_1_widget = widgets.Dropdown(
            options=comparisons_options,
            description="Comparison 1:",
            value=comparisons_options[0],
            layout=widgets.Layout(width="30%"),
            style={"description_width": "100px"},
        )
        # Create a widget to define the selection: <=, >=, <, >, abs() >=, abs() <=
        orientation_1_widget = widgets.Dropdown(
            options=["<=", ">=", "<", ">", "abs() >=", "abs() <="],
            description="",
            value=">=",
            layout=widgets.Layout(width="10%"),
        )
        # Create a widget to define the value of the selection
        value_1_widget = widgets.FloatText(
            value=0,
            description="",
            layout=widgets.Layout(width="10%"),
        )
        # Create a color widget for the condition 1
        color_1_widget = widgets.ColorPicker(
            concise=True,
            description="",
            value="blue",
            # layout=widgets.Layout(width="10%"),
        )
        # Create a checkbox to let user decide if they want to highlight genes with a specific value for the condition 1
        highlight_1_widget = widgets.Checkbox(
            value=True,
            description="Highlight genes?",
        )

        # Create a condition 2 widget
        comparison_2_widget = widgets.Dropdown(
            options=comparisons_options,
            description="Comparison 2:",
            value=comparisons_options[0],
            layout=widgets.Layout(width="30%"),
            style={"description_width": "100px"},
        )
        # Create a widget to define the selection: <=, >=, <, >, abs() >=, abs() <=
        orientation_2_widget = widgets.Dropdown(
            options=["<=", ">=", "<", ">", "abs() >=", "abs() <="],
            description="",
            value=">=",
            layout=widgets.Layout(width="10%"),
        )
        # Create a widget to define the value of the selection
        value_2_widget = widgets.FloatText(
            value=0,
            description="",
            layout=widgets.Layout(width="10%"),
        )
        # Create a color widget for the condition 2
        color_2_widget = widgets.ColorPicker(concise=True, description="", value="red")
        # Create a checkbox to let user decide if they want to highlight genes with a specific value for the condition 1
        highlight_2_widget = widgets.Checkbox(
            value=True,
            description="Highlight genes?",
        )

        # Create a color widget for elements identified for both conditions
        color_intersection_widget = widgets.ColorPicker(
            concise=True, description="Intersection", value="purple"
        )

        # Create a tool widget
        tools_widget = widgets.Dropdown(
            options=tools_options,
            description="Tool:",
            style={"description_width": "100px"},
        )

        return (
            tools_widget,
            widgets.HBox(
                [
                    comparison_1_widget,
                    orientation_1_widget,
                    value_1_widget,
                    color_1_widget,
                    highlight_1_widget,
                ]
            ),
            widgets.HBox(
                [
                    comparison_2_widget,
                    orientation_2_widget,
                    value_2_widget,
                    color_2_widget,
                    highlight_2_widget,
                ]
            ),
            color_intersection_widget,
            widgets.HBox([gene_selection_widget, color_gene_widget]),
        )

    parameters_widgets = widgets.VBox(parameters_widgets())

    # Display the widgets
    display(parameters_widgets)

    # Create a button to plot the data
    plot_button = widgets.Button(description="Plot")
    display(plot_button)

    # On button click, plot the data
    @plot_button.on_click
    def plot_button_clicked(b):
        """Plot the data."""

        def plot_comparison(data, column_1, column_2, elements_column, color_dict):
            """Plot the comparison using altair."""

            # Définir les limites de vos données
            min_value = min(data[column_1].min(), data[column_2].min())
            max_value = max(data[column_1].max(), data[column_2].max())

            # Créer une échelle avec les mêmes limites pour les axes x et y
            scale = alt.Scale(domain=(min_value, max_value))

            range_color = list(color_dict.values())
            domain_color = list(color_dict.keys())

            # Add none to the range color
            range_color.append("grey")
            domain_color.append("Others")

            # Superposer les lignes sur le graphique
            chart = (
                (
                    (
                        alt.Chart(data)
                        .transform_calculate(
                            order="datum.color == 'Others' ? 0 : (datum.color == 'Selection' ? 2 : 1)"
                        )
                        .mark_circle()
                        .encode(
                            x=alt.X(
                                column_1,
                                scale=scale,
                                title=column_1,
                                axis=alt.Axis(
                                    values=list(
                                        np.arange(min_value, max_value + 0.5, 0.5)
                                    )
                                ),
                            ),
                            y=alt.Y(
                                column_2,
                                scale=scale,
                                title=column_2,
                                axis=alt.Axis(
                                    values=list(
                                        np.arange(min_value, max_value + 0.5, 0.5)
                                    )
                                ),
                            ),
                            tooltip=[elements_column, column_1, column_2],
                            # Color the points in blue if they are selected
                            color=alt.Color(
                                "color:N",
                                scale=alt.Scale(
                                    range=range_color,
                                    domain=domain_color,
                                ),
                                sort="ascending",
                                # Set the legend title
                                legend=alt.Legend(title="Highlighted elements:"),
                            ),
                            order="order:O",
                            # Set opacity to 1.0 for all points
                            opacity=alt.value(1.0),
                        )
                        .interactive()
                    )
                    # Add diagonal line
                    + alt.Chart(pd.DataFrame({"x": [min_value * 2, max_value * 2]}))
                    .mark_line(color="#636363")
                    .encode(x="x", y="x")
                    # Add a vertical line to highlight the x = 0 line, in black
                    + alt.Chart(pd.DataFrame({"x": [0, 0]}))
                    .mark_rule(color="#636363")
                    .encode(x="x")
                    # Add a horizontal line to highlight the y = 0 line, in black
                    + alt.Chart(pd.DataFrame({"y": [0, 0]}))
                    .mark_rule(color="#636363")
                    .encode(y="y")
                    + alt.Chart(data.query("selected == True"))
                    .mark_text(dy=10, dx=20, color=selected_genes_color)
                    .encode(
                        x=column_1,
                        y=column_2,
                        text=elements_column,
                    )
                )
                .configure_axis(grid=False)
                .configure_view(stroke=None)
                .properties(width=500, height=500)
            )

            if "abs" in parameters_widgets.children[1].children[1].value:
                chart = chart + alt.Chart(
                    pd.DataFrame(
                        {
                            "x": [
                                parameters_widgets.children[1].children[2].value,
                                -parameters_widgets.children[1].children[2].value,
                            ]
                        }
                    )
                ).mark_rule(color="grey", strokeDash=[3, 3]).encode(x="x")
            else:
                chart = chart + alt.Chart(
                    pd.DataFrame(
                        {"x": [parameters_widgets.children[1].children[2].value]}
                    )
                ).mark_rule(color="grey", strokeDash=[3, 3]).encode(x="x")
            if "abs" in parameters_widgets.children[2].children[1].value:
                chart = chart + alt.Chart(
                    pd.DataFrame(
                        {
                            "y": [
                                parameters_widgets.children[2].children[2].value,
                                -parameters_widgets.children[2].children[2].value,
                            ]
                        }
                    )
                ).mark_rule(color="grey", strokeDash=[3, 3]).encode(y="y")
            else:
                chart = chart + alt.Chart(
                    pd.DataFrame(
                        {"y": [parameters_widgets.children[2].children[2].value]}
                    )
                ).mark_rule(color="grey", strokeDash=[3, 3]).encode(y="y")

            display(chart)

        def _get_elements_name_column_by_tool(tool):
            """Get the column name for the elements by tool."""
            columns_by_tool = {
                "MAGeCK_MLE": "Gene",
                "MAGeCK_RRA": "id",
                "CRISPhieRmix": "gene",
                "SSREA": "pathway",
                "directional_scoring_method": "Gene",
            }
            return columns_by_tool[tool]

        def _get_score_columns_by_tool(tool, condition):
            """Get the columns to plot by tool."""
            columns_by_tool = {
                "MAGeCK_MLE": f"{condition}|beta",
                "MAGeCK_RRA": "neg|lfc",
                "CRISPhieRmix": "mean_log2FoldChange",
                "SSREA": "NES",
                "directional_scoring_method": "score",
            }
            return columns_by_tool[tool]

        # with output:
        # Clear the output
        # output.clear_output()
        # Get the parameters
        tool = parameters_widgets.children[0].value

        # Comparison 1 parameters
        comparison_1 = parameters_widgets.children[1].children[0].value
        orientation_1 = parameters_widgets.children[1].children[1].value
        value_1 = parameters_widgets.children[1].children[2].value
        color_1 = parameters_widgets.children[1].children[3].value
        highlight_1 = parameters_widgets.children[1].children[4].value

        # Comparison 2 parameters
        comparison_2 = parameters_widgets.children[2].children[0].value
        orientation_2 = parameters_widgets.children[2].children[1].value
        value_2 = parameters_widgets.children[2].children[2].value
        color_2 = parameters_widgets.children[2].children[3].value
        highlight_2 = parameters_widgets.children[2].children[4].value

        # Get the color for the intersection
        color_intersection = parameters_widgets.children[3].value

        # Get the data
        comparison_1_data, comparison_2_data = get_data(
            tool, comparison_1, comparison_2, results_directory
        )

        # If both comparisons are the same, display an error message
        if comparison_1 == comparison_2:
            display_warning("Please choose two different comparisons.")
            # print("Please choose two different comparisons.")
            return

        # Get the selected genes and their color
        selected_genes = parameters_widgets.children[4].children[0].value.split(",")
        selected_genes_color = parameters_widgets.children[4].children[1].value

        # Get the columns to plot
        column_1 = _get_score_columns_by_tool(tool, comparison_1.split("_vs_")[0])
        column_2 = _get_score_columns_by_tool(tool, comparison_2.split("_vs_")[0])
        # Get column with elements name
        elements_column = _get_elements_name_column_by_tool(tool)

        # Combine the data: name, comparison_1, comparison_2
        combined_data = comparison_1_data[[elements_column, column_1]].merge(
            comparison_2_data[[elements_column, column_2]],
            on=elements_column,
            suffixes=("_" + comparison_1, "_" + comparison_2),
            # Force suffixes to be added to the columns
            how="outer",
        )

        if not column_1 + "_" + comparison_1 in combined_data.columns:
            # Add a suffix to the column name
            combined_data[column_1 + "_" + comparison_1] = comparison_1_data[column_1]
        column_1 = column_1 + "_" + comparison_1
        if not column_2 + "_" + comparison_2 in combined_data.columns:
            # Add a suffix to the column name
            combined_data[column_2 + "_" + comparison_2] = comparison_2_data[column_2]
        column_2 = column_2 + "_" + comparison_2

        # Add an annotation column to distinguish the elements selected by the user
        combined_data["selected"] = combined_data[elements_column].isin(selected_genes)

        def mask_condition(data, column, orientation, value):
            """Mask the data based on the condition."""
            if orientation == "<=":
                mask = data[column] <= value
            elif orientation == ">=":
                mask = data[column] >= value
            elif orientation == "<":
                mask = data[column] < value
            elif orientation == ">":
                mask = data[column] > value
            elif orientation == "abs() <=":
                mask = abs(data[column]) <= value
            elif orientation == "abs() >=":
                mask = abs(data[column]) >= value
            return mask

        # If user selected to highlight genes, add a column to the dataframe for the genes passing the threshold
        combined_data["highlight_1"] = mask_condition(
            combined_data, column_1, orientation_1, value_1
        )
        combined_data["highlight_2"] = mask_condition(
            combined_data, column_2, orientation_2, value_2
        )

        conditions = [
            combined_data["selected"],
            combined_data["highlight_1"] & combined_data["highlight_2"],
            combined_data["highlight_1"],
            combined_data["highlight_2"],
        ]

        label_1 = f"{comparison_1} {orientation_1} {value_1}"
        label_2 = f"{comparison_2} {orientation_2} {value_2}"

        outputs = [
            "Selection",
            "Intersection",
            label_1,
            label_2,
        ]

        # If no highlight is selected, remove the first element of the lists
        to_remove = []
        if not highlight_1 or not highlight_2:
            to_remove.append(1)
        if not highlight_1:
            to_remove.append(2)
        if not highlight_2:
            to_remove.append(3)

        conditions = [
            conditions[i] for i in range(len(conditions)) if i not in to_remove
        ]
        outputs = [outputs[i] for i in range(len(outputs)) if i not in to_remove]

        # Create a 'color' column.
        # If highlight_1 is True, set the color to color_1
        # If highlight_2 is True, set the color to color_2
        # If both are True, set the color to color_intersection
        # If selected is True, set the color to selected_genes_color
        # Else, set the color to grey
        combined_data["color"] = np.select(
            conditions,
            outputs,
            default="Others",
        )

        # display(combined_data)

        # display(combined_data.groupby("color").size())

        color_dict = {
            # Highlight the genes passing the threshold 1
            label_1: color_1,
            # Highlight the genes passing the threshold 2
            label_2: color_2,
            # Highlight the genes in the intersection
            "Intersection": color_intersection,
            # Genes selected by the user
            "Selection": selected_genes_color,
            # Others genes
            "Others": "grey",
        }

        # display(combined_data)

        # Plot the data
        plot_comparison(combined_data, column_1, column_2, elements_column, color_dict)

        # Create a dropdown to select the gene sets
        enrichr_bases = get_enrichr_bases()
        enrichr_bases_widget = widgets.SelectMultiple(
            options=enrichr_bases,
            description="EnrichR genesets:",
            rows=12,
            layout=widgets.Layout(width="400px"),
            style={"description_width": "initial"},
        )

        display(enrichr_bases_widget)

        label_to_column = {
            label_1: "highlight_1",
            label_2: "highlight_2",
        }

        # # Retrieve genes from each 'color' category and create a text area for each category
        # for category in color_dict:
        #     if category not in ["Selection", "Others"]:
        #     genes = combined_data[combined_data["color"] == category][elements_column]

        for category in [label_1, label_2, "Intersection"]:
            if category == "Intersection":
                genes = combined_data[
                    (combined_data["highlight_1"]) & (combined_data["highlight_2"])
                ][elements_column]
            elif category in [label_1, label_2]:
                genes = combined_data[combined_data[label_to_column[category]]][
                    elements_column
                ]
            textarea = widgets.Textarea(
                value="\n".join(genes),
                description=f"{category} genes:",
                disabled=True,
                layout=widgets.Layout(height="200px", width="auto"),
                style={"description_width": "200px"},
            )

            # Create a button to run enrichr on the genes
            enrichr_button = widgets.Button(
                description=f"Run EnrichR on {category}",
                # Increase button width to match description width
                layout=widgets.Layout(width="auto"),
            )
            display(widgets.HBox([textarea, enrichr_button]))

            def enrichr_button_clicked(b):
                """Run enrichr on the genes."""
                category = b.description[15:]
                if category == "Intersection":
                    genes = combined_data[
                        (combined_data["highlight_1"]) & (combined_data["highlight_2"])
                    ][elements_column]
                elif category in [label_1, label_2]:
                    genes = combined_data[combined_data[label_to_column[category]]][
                        elements_column
                    ]
                if len(genes) == 0:
                    display_warning(f"No genes in the category {category}.")
                    # print(f"No genes in the category {category}.")
                    return
                print(f"Running EnrichR on <{category}> genes.")
                print(f"Number of genes: {len(genes)}")
                if not enrichr_bases_widget.value:
                    display_warning("Please select at least one geneset above.")
                    # print("Please select at least one gene set.")
                    return
                for base in enrichr_bases_widget.value:
                    enrichr_res = get_enrichr_results(
                        genes,
                        category,
                        base,
                    )
                    table = create_enrichr_table(enrichr_res)
                    chart = enrichmentBarPlot(
                        table,
                        10,
                        f"EnrichR results for {category} genes",
                        "red",
                        "blue",
                        base,
                    )
                    display(chart)

            enrichr_button.on_click(enrichr_button_clicked)

        # Display the text in black
        display(
            HTML(
                """
            <style>
                .widget-textarea textarea:disabled {
                    color: black !important;
                    opacity: 1 !important;
                }
            </style>
            """
            )
        )

        # # Display the text areas in a Hbox
        # display(widgets.HBox([upper_left_textarea, upper_right_textarea]))
        # display(widgets.HBox([lower_left_textarea, lower_right_textarea]))
