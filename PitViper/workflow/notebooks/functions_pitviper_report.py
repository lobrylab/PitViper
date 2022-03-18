import json
import os
import re
import warnings
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
from ipywidgets import interact
from rpy2.rinterface import RRuntimeWarning
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.lib.dplyr import DataFrame
from rpy2.robjects.packages import importr
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import zscore
from sklearn import datasets, decomposition

warnings.filterwarnings("ignore")

warnings.filterwarnings("ignore", category=RRuntimeWarning)

depmap = importr("depmap")
experimentHub = importr("ExperimentHub")
utils = importr("utils")

## Remove Altair max rows
alt.data_transformers.disable_max_rows()

pd.options.mode.chained_assignment = None


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
    print("Working directory: ", os.getcwd())

    with open(output, "w") as out:
        print("Notebook was runned.", file=out)


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


def import_results(token: str):
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
    if content["mageck_rra_activate"] == "True":
        tools.append("MAGeCK_RRA")
    if content["bagel_activate"] == "True":
        tools.append("BAGEL")
    if content["crisphiermix_activate"] == "True":
        tools.append("CRISPhieRmix")
    if content["filtering_activate"] == "True":
        tools.append("in_house_method")
    if content["gsea_activate"] == "True":
        tools.append("GSEA-like")

    results_directory = "results/%s/" % token
    print("Results directory: %s \n" % results_directory)

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
                if file.endswith(".txt") or file.endswith(".bf"):
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
    if "GSEA-like" in tools_available:
        for condition in tools_available["GSEA-like"]:
            treatment = condition.split("_vs_")[0]
            for file_suffixe in ["%s_all-elements_GSEA-like.txt"]:
                table = tools_available["GSEA-like"][condition][
                    file_suffixe % condition
                ]
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
        print("No mapping QC file to show.")
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
    return s


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
    path_qc = content["normalized_count_table"]
    if not path.exists(path_qc):
        print("No count file to show.")
        return 0
    table = pd.read_csv(path_qc, sep="\t")

    table.iloc[:, 2:] = table.iloc[:, 2:] + 1
    table.iloc[:, 2:] = table.iloc[:, 2:].apply(np.log2)

    chart = (
        alt.Chart(table)
        .transform_fold(list(table.columns[2:]), as_=["Measurement_type", "counts"])
        .transform_density(
            density="counts",
            bandwidth=0.3,
            groupby=["Measurement_type"],
            extent=[0, 20],
            counts=True,
            steps=200,
        )
        .mark_line()
        .encode(
            alt.X("value:Q", axis=alt.Axis(title="log2(read count)")),
            alt.Y("density:Q"),
            alt.Color("Measurement_type:N"),
            tooltip=["Measurement_type:N", "value:Q", "density:Q"],
        )
        .properties(width=width, height=height)
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

    TSV = pd.read_csv(content["tsv_file"], sep="\t")
    cts_file = content["normalized_count_table"]
    cts = pd.read_csv(cts_file, sep="\t")
    X = cts[cts.columns[2:]].to_numpy().T
    d = dict(zip(TSV.replicate, TSV.condition))
    y = [d[k] for k in cts.columns[2:]]
    y = np.array(y)
    y_bis = np.array(cts.columns[2:])

    pca = decomposition.PCA(n_components=2)
    pca.fit(X)
    X = pca.transform(X)

    a = pd.DataFrame(X, columns=["PC1", "PC2"])
    b = pd.DataFrame(y, columns=["condition"])
    c = pd.DataFrame(y_bis, columns=["replicate"])

    df_c = pd.concat([a, b, c], axis=1)

    source = df_c

    PC1_explained_variance_ratio = round(pca.explained_variance_ratio_[0] * 100, 2)
    PC2_explained_variance_ratio = round(pca.explained_variance_ratio_[1] * 100, 2)

    pca_2d = (
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

    return pca_2d


def getEnrichrResults(genes: list, description: str, gene_set_library: list):
    """Get and return EnrichR results for a given list of genes.

    Args:
        genes (list): List of genes to use for EnrichR analysis.
        description (str): Description of the genes list.
        gene_set_library (list): List of genes set to use.

    Raises:
        Exception: Error analyzing gene list.
        Exception: Error fetching enrichment results.

    Returns:
        dict: EnrichR results
    """
    ENRICHR_URL = "http://maayanlab.cloud/Enrichr/addList"
    genes_str = "\n".join(genes)
    description = description
    payload = {"list": (None, genes_str), "description": (None, description)}

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception("Error analyzing gene list")

    data = json.loads(response.text)

    userListId = data["userListId"]

    ENRICHR_URL = "http://maayanlab.cloud/Enrichr/enrich"
    query_string = "?userListId=%s&backgroundType=%s"
    user_list_id = userListId
    gene_set_library = gene_set_library[:-1]
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
    )
    if not response.ok:
        raise Exception("Error fetching enrichment results")

    data = json.loads(response.text)
    return data


def createEnrichrTable(enrichrResults: dict):
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
            title=description + " (%s)" % base[:-1],
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


def GSEA_like_data(
    comparison="", control="", tool="", results_directory="", tools_available=""
):
    """Return GSEA-like results as pandas dataframe."""
    tables_list = []
    for _comparison in os.listdir(os.path.join(results_directory, tool)):
        if _comparison.split("_vs_")[-1] == control:
            keys_list = list(tools_available[tool][_comparison].keys())
            for key in keys_list:
                if key.endswith("_all-elements_GSEA-like.txt"):
                    break
            data = tools_available[tool][_comparison][key]
            trt = _comparison.split("_vs_")[0]
            data["condition"] = trt
            tables_list.append(data)
        if _comparison == comparison:
            data = tools_available[tool][_comparison][
                "%s_all-elements_GSEA-like.txt" % comparison
            ]
            tables_list.append(data)
    result = pd.concat(tables_list)
    return result


def in_house_method_data(
    comparison="", control="", tool="", results_directory="", tools_available=""
):
    """Return in-house method results as pandas dataframe."""
    tables_list = []
    for _comparison in os.listdir(os.path.join(results_directory, tool)):
        if _comparison.split("_vs_")[-1] == control:
            keys_list = list(tools_available[tool][_comparison].keys())
            for key in keys_list:
                if key.endswith("_all-elements_in-house.txt"):
                    break
            data = tools_available[tool][_comparison][key]
            trt = _comparison.split("_vs_")[0]
            data["condition"] = trt
            tables_list.append(data)
        if _comparison == comparison:
            data = tools_available[tool][_comparison][
                "%s_all-elements_in-house.txt" % comparison
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
        print("ERROR.")
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
    comparison="", control="", tool="BAGEL", results_directory="", tools_available=""
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

    tables_list = []
    if control != "" and comparison == "":
        mode = True
    elif comparison != "" and control == "":
        mode = False
    else:
        print("ERROR.")
    for _comparison in os.listdir(os.path.join(results_directory, tool)):
        if _comparison.split("_vs_")[-1] == control:
            keys_list = list(tools_available[tool][_comparison].keys())
            for key in keys_list:
                if key.endswith("_BAGEL_output.bf"):
                    break
            data = tools_available[tool][_comparison][key]
            trt = _comparison.split("_vs_")[0]
            data["condition"] = trt
            tables_list.append(data)
        if _comparison == comparison:
            data = tools_available[tool][_comparison][
                "%s_BAGEL_output.bf" % _comparison
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


def tool_results(results_directory, tools_available, token):
    """Display selected method's results for all genes."""
    config = "./config/%s.yaml" % token
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
        min=0.0, max=1.0, step=0.01, value=0.05, description="FDR cut-off"
    )

    color_sig_widget = widgets.ColorPicker(
        concise=False, description="Significant color:", value="red"
    )
    color_non_widget = widgets.ColorPicker(
        concise=False, description="Non-significant color:", value="gray"
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
        tool = "MAGeCK_MLE"
        significant_label = "FDR < %s" % fdr_cutoff
        non_significant_label = "FDR ≥ %s" % fdr_cutoff
        highlight_label = "Gene(s) of interest"
        treatment, control = comparison.split("_vs_")
        source = MAGeCK_MLE_data(
            comparison=comparison,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )
        source["default_rank"] = source[treatment + "|beta"].rank()
        source.loc[
            source[treatment + "|fdr"] < fdr_cutoff, "significant"
        ] = significant_label
        source.loc[
            source[treatment + "|fdr"] >= fdr_cutoff, "significant"
        ] = non_significant_label
        source.loc[source.Gene.isin(elements), "significant"] = highlight_label
        domain = [significant_label, non_significant_label, highlight_label]
        range_ = [sig, non_sig, "blue"]
        chart = (
            alt.Chart(source, title="MAGeCK MLE (%s)" % comparison)
            .mark_circle(size=60)
            .encode(
                x=alt.X("default_rank:Q", axis=alt.Axis(title="Rank")),
                y=alt.Y(
                    treatment + "|beta:Q", axis=alt.Axis(title="%s beta" % treatment)
                ),
                tooltip=[
                    "Gene",
                    "sgRNA",
                    treatment + "|beta",
                    treatment + "|fdr",
                    "significant",
                    "default_rank",
                ],
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                order=alt.Order("significant:N"),
            )
            .properties(width=800, height=400)
            .interactive()
        )
        line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule().encode(y="y")
        text = (
            alt.Chart(source.query("significant == 'Gene(s) of interest'"))
            .mark_text(dy=10, dx=20, color="blue")
            .encode(
                x=alt.X("default_rank:Q"),
                y=alt.Y(treatment + "|beta:Q"),
                text=alt.Text("Gene"),
            )
        )
        chart = chart + line + text
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
        tool = "MAGeCK_RRA"
        significant_label = "FDR < %s" % fdr_cutoff
        non_significant_label = "FDR ≥ %s" % fdr_cutoff
        highlight_label = "Gene(s) of interest"
        treatment, control = comparison.split("_vs_")
        source = MAGeCK_RRA_data(
            comparison=comparison,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )
        source["default_rank"] = source["neg|lfc"].rank(ascending=True)
        source.loc[source["neg|fdr"] < fdr_cutoff, "significant"] = significant_label
        source.loc[
            source["neg|fdr"] >= fdr_cutoff, "significant"
        ] = non_significant_label
        source.loc[source.id.isin(elements), "significant"] = highlight_label
        domain = [significant_label, non_significant_label, highlight_label]
        range_ = [sig, non_sig, "blue"]
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
                    "neg|fdr",
                    "significant",
                    "neg|rank",
                    "default_rank",
                ],
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                order=alt.Order("significant:N"),
            )
            .properties(width=800, height=400)
            .interactive()
        )
        line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule().encode(y="y")
        text = (
            alt.Chart(source.query("significant == 'Gene(s) of interest'"))
            .mark_text(dy=10, dx=20, color="blue")
            .encode(
                x=alt.X("default_rank:Q"), y=alt.Y("neg|lfc:Q"), text=alt.Text("id")
            )
        )

        chart = chart + line + text
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
        tool = "CRISPhieRmix"
        significant_label = "FDR < %s" % fdr_cutoff
        non_significant_label = "FDR ≥ %s" % fdr_cutoff
        highlight_label = "Gene(s) of interest"
        treatment, control = comparison.split("_vs_")
        source = CRISPhieRmix_data(
            comparison=comparison,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )
        source["default_rank"] = source["top_3_mean_log2FoldChange"].rank(
            method="dense"
        )
        source.loc[source["locfdr"] < fdr_cutoff, "significant"] = significant_label
        source.loc[
            source["locfdr"] >= fdr_cutoff, "significant"
        ] = non_significant_label
        source.loc[source.gene.isin(elements), "significant"] = highlight_label
        domain = [significant_label, non_significant_label, highlight_label]
        range_ = [sig, non_sig, "blue"]
        chart = (
            alt.Chart(source, title="CRISPhieRmix (%s)" % comparison)
            .mark_circle(size=60)
            .encode(
                x=alt.X("default_rank:Q", axis=alt.Axis(title="Rank")),
                y=alt.Y(
                    "top_3_mean_log2FoldChange:Q",
                    axis=alt.Axis(
                        title="%s top-3 sgRNAs log2FoldChange average" % treatment
                    ),
                ),
                tooltip=[
                    "gene",
                    "locfdr",
                    "FDR",
                    "significant",
                    "default_rank",
                    "top_3_mean_log2FoldChange",
                ],
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                order=alt.Order("significant:N"),
            )
            .properties(width=800, height=400)
            .interactive()
        )
        line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule().encode(y="y")
        text = (
            alt.Chart(source.query("significant == 'Gene(s) of interest'"))
            .mark_text(dy=10, dx=20, color="blue")
            .encode(
                x=alt.X("default_rank:Q"),
                y=alt.Y("top_3_mean_log2FoldChange:Q"),
                text=alt.Text("gene"),
            )
        )
        chart = chart + line + text
        display(chart)

    def _in_house_snake_plot(
        comparison,
        fdr_cutoff,
        non_sig,
        sig,
        results_directory,
        tools_available,
        elements,
    ):
        tool = "in_house_method"
        significant_label = "Filter: pass"
        non_significant_label = "Filter: don't pass"
        highlight_label = "Gene(s) of interest"
        treatment, control = comparison.split("_vs_")
        source = in_house_method_data(
            comparison=comparison,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )
        source["default_rank"] = source["score"].rank(method="first")
        source.loc[
            source.category.isin(["down", "up"]), "significant"
        ] = significant_label
        source.loc[
            ~source.category.isin(["down", "up"]), "significant"
        ] = non_significant_label
        source.loc[source.Gene.isin(elements), "significant"] = highlight_label
        domain = [significant_label, non_significant_label, highlight_label]
        range_ = [sig, non_sig, "blue"]
        chart = (
            alt.Chart(source, title="In-house method (%s)" % comparison)
            .mark_circle(size=60)
            .encode(
                x=alt.X("default_rank:Q", axis=alt.Axis(title="Rank")),
                y=alt.Y("score:Q", axis=alt.Axis(title="%s score" % treatment)),
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
                order=alt.Order("significant:N"),
            )
            .properties(width=800, height=400)
            .interactive()
        )
        line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule().encode(y="y")
        text = (
            alt.Chart(source.query("significant == 'Gene(s) of interest'"))
            .mark_text(dy=10, dx=20, color="blue")
            .encode(
                x=alt.X("default_rank:Q"), y=alt.Y("score:Q"), text=alt.Text("Gene")
            )
        )
        chart = chart + line + text
        display(chart)

    def _GSEA_like_snake_plot(
        comparison,
        fdr_cutoff,
        non_sig,
        sig,
        results_directory,
        tools_available,
        elements,
    ):
        tool = "GSEA-like"
        significant_label = "FDR < %s" % fdr_cutoff
        non_significant_label = "FDR ≥ %s" % fdr_cutoff
        highlight_label = "Gene(s) of interest"
        treatment, control = comparison.split("_vs_")
        source = GSEA_like_data(
            comparison=comparison,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )
        source["default_rank"] = source[["NES"]].rank(method="dense")
        source.loc[abs(source["padj"]) < fdr_cutoff, "significant"] = significant_label
        source.loc[
            abs(source["padj"]) >= fdr_cutoff, "significant"
        ] = non_significant_label
        source.loc[source.pathway.isin(elements), "significant"] = highlight_label
        domain = [significant_label, non_significant_label, highlight_label]
        range_ = [sig, non_sig, "blue"]

        chart = (
            alt.Chart(source, title="GSEA-like method (%s)" % comparison)
            .mark_circle(size=60)
            .encode(
                x=alt.X("default_rank:Q", axis=alt.Axis(title="Rank")),
                y=alt.Y("NES:Q", axis=alt.Axis(title="%s NES" % treatment)),
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
                order=alt.Order("significant:N"),
            )
            .properties(width=800, height=400)
            .interactive()
        )
        line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule().encode(y="y")
        text = (
            alt.Chart(source.query("significant == 'Gene(s) of interest'"))
            .mark_text(dy=10, dx=20, color="blue")
            .encode(
                x=alt.X("default_rank:Q"), y=alt.Y("NES:Q"), text=alt.Text("pathway")
            )
        )
        chart = chart + line + text
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
        tool = "BAGEL"
        significant_label = "BF > 0"
        non_significant_label = "BF ≤ 0"
        highlight_label = "Gene(s) of interest"
        treatment, control = comparison.split("_vs_")
        source = BAGEL_data(
            comparison=comparison,
            control="",
            tool=tool,
            results_directory=results_directory,
            tools_available=tools_available,
        )
        source["default_rank"] = source["BF"].rank(method="dense", ascending=False)
        source.loc[source["BF"] > 0, "significant"] = significant_label
        source.loc[source["BF"] <= 0, "significant"] = non_significant_label
        source.loc[source.GENE.isin(elements), "significant"] = highlight_label
        domain = [significant_label, non_significant_label, highlight_label]
        range_ = [sig, non_sig, "blue"]

        chart = (
            alt.Chart(source, title="BAGEL (%s)" % comparison)
            .mark_circle(size=60)
            .encode(
                x=alt.X("default_rank:Q", axis=alt.Axis(title="Rank")),
                y=alt.Y("BF:Q", axis=alt.Axis(title="%s Bayesian Factor" % treatment)),
                tooltip=["GENE", "BF", "STD", "significant", "default_rank"],
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                order=alt.Order("significant:N"),
            )
            .properties(width=800, height=400)
            .interactive()
        )
        line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule().encode(y="y")
        text = (
            alt.Chart(source.query("significant == 'Gene(s) of interest'"))
            .mark_text(dy=10, dx=20, color="blue")
            .encode(x=alt.X("default_rank:Q"), y=alt.Y("BF:Q"), text=alt.Text("GENE"))
        )
        chart = chart + line + text
        display(chart)

    def _plot(event):
        tool = tools_widget.value
        comparison = comparisons_widget.value
        fdr_cutoff = fdr_widget.value
        non_sig = color_non_widget.value
        sig = color_sig_widget.value
        elements = element.value.split(",")
        if "MAGeCK_RRA" in tool:
            _MAGeCK_RRA_snake_plot(
                comparison,
                fdr_cutoff,
                non_sig,
                sig,
                results_directory,
                tools_available,
                elements,
            )
        if "MAGeCK_MLE" in tool:
            _MAGeCK_MLE_snake_plot(
                comparison,
                fdr_cutoff,
                non_sig,
                sig,
                results_directory,
                tools_available,
                elements,
            )
        if "CRISPhieRmix" in tool:
            _CRISPhieRmix_snake_plot(
                comparison,
                fdr_cutoff,
                non_sig,
                sig,
                results_directory,
                tools_available,
                elements,
            )
        if "in_house_method" in tool:
            _in_house_snake_plot(
                comparison,
                fdr_cutoff,
                non_sig,
                sig,
                results_directory,
                tools_available,
                elements,
            )
        if "GSEA-like" in tool:
            _GSEA_like_snake_plot(
                comparison,
                fdr_cutoff,
                non_sig,
                sig,
                results_directory,
                tools_available,
                elements,
            )
        if "BAGEL" in tool:
            _BAGEL_snake_plot(
                comparison,
                fdr_cutoff,
                non_sig,
                sig,
                results_directory,
                tools_available,
                elements,
            )
        else:
            print("Choose a tool.")

    display(tools_widget)
    display(element)
    display(comparisons_widget)
    display(fdr_widget)
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
                print("%s was not found in %s." % (gene, cts_file))
            else:
                gene_cts = cts.loc[cts.Gene == gene]
                gene_cts = gene_cts.loc[gene_cts.variable.isin(conditions.value)]
                z_scores = gene_cts.groupby(["sgRNA"]).value.transform(
                    lambda x: zscore(x, ddof=1)
                )
                gene_cts["z-score"] = z_scores
                sum_by_group = (
                    gene_cts.groupby(["sgRNA"]).agg({"value": "sum"}).reset_index()
                )
                sgRNA_to_keep = sum_by_group[sum_by_group.value > 50].sgRNA.values
                gene_cts = gene_cts[gene_cts.sgRNA.isin(sgRNA_to_keep)]
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

    display(element)
    display(conditions)

    button = widgets.Button(description="Show!")

    display(button)

    def show_plot(source, sort_cols, gene):
        selection = alt.selection_multi(fields=["sgRNA"], bind="legend")
        base = (
            alt.Chart(source, title="%s replicates normalized read counts" % gene)
            .transform_fold(sort_cols)
            .encode(
                x=alt.X("key:N", sort=sort_cols, axis=alt.Axis(title="Condition")),
                y=alt.Y("value:Q", axis=alt.Axis(title="Condition")),
                color="sgRNA:N",
                detail="replicate_group:N",
                tooltip=["sgRNA", "value:Q", "replicate_group"],
            )
        )
        points = (
            base.mark_circle()
            .encode(opacity=alt.condition(selection, alt.value(1), alt.value(0.0)))
            .add_selection(selection)
            .properties(width=600)
        )
        lines = base.mark_line().encode(
            opacity=alt.condition(selection, alt.value(1), alt.value(0.0))
        )
        chart = points + lines
        display(chart.interactive())

    def on_button_clicked(b):
        # Read normalized counts file
        cts = pd.read_csv(cts_file, sep="\t")
        # Transform dataframe: one value per line
        cts = pd.melt(cts, id_vars=["sgRNA", "Gene"])
        genes = element.value.split(",")
        for gene in genes:
            if not gene in list(cts.Gene):
                print("Element '%s' was not found in %s." % (gene, cts_file))
            else:
                sort_cols = conditions.value
                gene_cts = cts.loc[cts.Gene == gene]
                source = gene_cts
                source = pd.merge(
                    source, design, left_on="variable", right_on="replicate"
                )
                boolean_series = source.condition.isin(sort_cols)
                source = source[boolean_series]
                source["replicate_group"] = source["replicate"].str.extract(
                    r"([1-9]+)$"
                )
                source = pd.pivot_table(
                    source,
                    values="value",
                    index=["sgRNA", "replicate_group"],
                    columns=["condition"],
                )
                source = source.reset_index()
                show_plot(source, sort_cols, gene)

    button.on_click(on_button_clicked)


def tool_results_by_element(results_directory, tools_available, token):
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
        elif tool == "GSEA-like":
            result = GSEA_like_data(
                comparison="",
                control=control.value,
                tool=tool,
                results_directory=results_directory,
                tools_available=tools_available,
            )
        elif tool == "in_house_method":
            result = in_house_method_data(
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
            elements_list = list(set(result.GENE))
        elif tool == "GSEA-like":
            elements_list = list(set(result.pathway))
        elif tool == "in_house_method":
            elements_list = list(set(result.Gene))
        return elements_list

    def update_control(update):
        ctrs = get_controls(results_directory, tools_available)
        control.options = ctrs

    def update_genes_list(update):
        result = CRISPhieRmix_data(
            comparison="",
            control=control.value,
            tool="CRISPhieRmix",
            results_directory=results_directory,
            tools_available=tools_available,
        )
        elements_list = result.gene.to_list()
        gene.options = elements_list
        gene.value = elements_list[0]

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

    tools_list = [tool for tool in list(tools_available.keys()) if tool != "DESeq2"]

    ctrs = get_controls(results_directory, tools_available, tools_list[0])
    control = widgets.Dropdown(
        options=ctrs, value=ctrs[0], description="Control:", disabled=False
    )
    control.observe(update_genes_list, "value")

    conditions = widgets.TagsInput(
        value=conditions_list, allowed_tags=conditions_list, allow_duplicates=False
    )

    fdr_cutoff = widgets.FloatSlider(
        description="FDR:", min=0.0, max=1.0, step=0.01, value=0.05
    )

    elements_list = get_genes_list(results_directory, tools_available, tools_list[0])
    gene = widgets.Combobox(
        placeholder="Choose one",
        options=elements_list,
        description="Element:",
        value=elements_list[0],
        ensure_option=False,
    )

    display(widgets.VBox([control, gene, fdr_cutoff, conditions]))

    def CRISPhieRmix_results(result, fdr_cutoff, control, gene, sort_cols):
        significant_label = "Yes"
        non_significant_label = "No"
        result.loc[result["locfdr"] < fdr_cutoff, "significant"] = significant_label
        result.loc[
            result["locfdr"] >= fdr_cutoff, "significant"
        ] = non_significant_label
        new_row = {
            "gene": gene,
            "condition": control,
            "significant": "Baseline",
            "locfdr": 1,
            "top_3_mean_log2FoldChange": 0,
        }
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.gene == gene]
        domain = [significant_label, non_significant_label, "Baseline"]
        range_ = ["red", "grey", "black"]
        plot = (
            alt.Chart(res)
            .transform_fold(sort_cols)
            .mark_circle(size=60)
            .mark_point(filled=True, size=100)
            .encode(
                y=alt.Y(
                    "top_3_mean_log2FoldChange",
                    axis=alt.Axis(title="Top-3 sgRNAs log2FoldChange average"),
                ),
                x=alt.X("condition:O", sort=sort_cols),
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                tooltip=["gene", "locfdr", "significant", "top_3_mean_log2FoldChange"],
            )
            .properties(title=gene + " (CRISPhieRmix)", width=100)
        )
        return plot

    def MAGeCK_RRA_results(result, fdr_cutoff, control, gene, sort_cols):
        significant_label = "Yes"
        non_significant_label = "No"
        result.loc[result["neg|fdr"] < fdr_cutoff, "significant"] = significant_label
        result.loc[
            result["neg|fdr"] >= fdr_cutoff, "significant"
        ] = non_significant_label
        new_row = {
            "id": gene,
            "condition": control,
            "significant": "Baseline",
            "neg|fdr": 1,
            "neg|lfc": 0,
        }
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.id == gene]
        domain = [significant_label, non_significant_label, "Baseline"]
        range_ = ["red", "grey", "black"]
        plot = (
            alt.Chart(res)
            .mark_circle(size=60)
            .mark_point(
                filled=True,
                size=100,
            )
            .encode(
                y=alt.Y("neg|lfc", axis=alt.Axis(title="FDR")),
                x=alt.X("condition:N", sort=sort_cols),
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                tooltip=["id", "neg|lfc", "significant", "neg|fdr", "condition"],
            )
            .properties(title=gene + " (MAGeCK RRA)", width=100)
        )
        return plot

    def MAGeCK_MLE_results(result, fdr_cutoff, control, gene, sort_cols):
        significant_label = "Yes"
        non_significant_label = "No"
        result.loc[result["fdr"] < fdr_cutoff, "significant"] = significant_label
        result.loc[result["fdr"] >= fdr_cutoff, "significant"] = non_significant_label
        new_row = {
            "Gene": gene,
            "condition": control,
            "significant": "Baseline",
            "fdr": 1,
            "beta": 0,
        }
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.Gene == gene]
        domain = [significant_label, non_significant_label, "Baseline"]
        range_ = ["red", "grey", "black"]
        plot = (
            alt.Chart(res)
            .mark_circle(size=60)
            .mark_point(
                filled=True,
                size=100,
            )
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

    def GSEA_like_results(result, fdr_cutoff, control, gene, sort_cols):
        significant_label = "Yes"
        non_significant_label = "No"
        result.loc[result["padj"] < fdr_cutoff, "significant"] = significant_label
        result.loc[result["padj"] >= fdr_cutoff, "significant"] = non_significant_label
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
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.pathway == gene]
        domain = [significant_label, non_significant_label, "Baseline"]
        range_ = ["red", "grey", "black"]
        plot = (
            alt.Chart(res)
            .mark_circle(size=60)
            .mark_point(
                filled=True,
                size=100,
            )
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
            .properties(title=gene + " (GSEA-like method)", width=100)
        )
        return plot

    def BAGEL_results(result, fdr_cutoff, control, gene, sort_cols):
        significant_label = "Yes"
        non_significant_label = "No"
        result.loc[result["BF"] > 0, "significant"] = significant_label
        result.loc[result["BF"] <= 0, "significant"] = non_significant_label
        new_row = {
            "GENE": gene,
            "condition": control,
            "significant": "Baseline",
            "BF": 0,
        }
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.GENE == gene]
        domain = [significant_label, non_significant_label, "Baseline"]
        range_ = ["red", "grey", "black"]
        plot = (
            alt.Chart(res)
            .mark_circle(size=60)
            .mark_point(
                filled=True,
                size=100,
            )
            .encode(
                y=alt.Y("BF", axis=alt.Axis(title="Bayesian factor")),
                x=alt.X("condition:N", sort=sort_cols),
                color=alt.Color(
                    "significant",
                    scale=alt.Scale(domain=domain, range=range_),
                    legend=alt.Legend(title="Significativity:"),
                ),
                tooltip=["GENE", "BF", "STD", "NumObs", "condition"],
            )
            .properties(title=gene + " (BAGEL)", width=100)
        )
        return plot

    def inhouse_method_results(result, fdr_cutoff, control, gene, sort_cols):
        significant_label = "Yes"
        non_significant_label = "No"
        result.loc[
            result.category.isin(["down", "up"]), "significant"
        ] = significant_label
        result.loc[
            ~result.category.isin(["down", "up"]), "significant"
        ] = non_significant_label
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
        domain = [significant_label, non_significant_label, "Baseline"]
        range_ = ["red", "grey", "black"]
        plot = (
            alt.Chart(res)
            .mark_circle(size=60)
            .mark_point(
                filled=True,
                size=100,
            )
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
            .properties(title=gene + " (Filtering method)", width=100)
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
                        fdr_cutoff.value,
                        control.value,
                        element,
                        conditions.value,
                    )
                elif tool == "MAGeCK_MLE":
                    plot = MAGeCK_MLE_results(
                        result,
                        fdr_cutoff.value,
                        control.value,
                        element,
                        conditions.value,
                    )
                elif tool == "GSEA-like":
                    plot = GSEA_like_results(
                        result,
                        fdr_cutoff.value,
                        control.value,
                        element,
                        conditions.value,
                    )
                elif tool == "in_house_method":
                    plot = inhouse_method_results(
                        result,
                        fdr_cutoff.value,
                        control.value,
                        element,
                        conditions.value,
                    )
                elif tool == "MAGeCK_RRA":
                    plot = MAGeCK_RRA_results(
                        result,
                        fdr_cutoff.value,
                        control.value,
                        element,
                        conditions.value,
                    )
                elif tool == "BAGEL":
                    plot = BAGEL_results(
                        result,
                        fdr_cutoff.value,
                        control.value,
                        element,
                        conditions.value,
                    )
                chart |= plot
            display(chart)

    button = widgets.Button(description="Show plot")
    display(button)
    button.on_click(on_button_clicked)


def enrichr_plots(token, pitviper_res):

    config = "./config/%s.yaml" % token
    content = open_yaml(config)
    if content["screen_type"] == "not_gene":
        return HTML(
            """<p style="color:red;background-color: white;padding: 0.5em;">This module is available only if genes symbol are available.</p>"""
        )
        # return "This module is available only if genes symbol are used."

    def update_conditions(update):
        conditions_list = list(pitviper_res[tool.value].keys())
        conditions.options = conditions_list
        conditions.value = conditions_list[0]

    BASES = open("workflow/notebooks/enrichr_list.txt", "r").readlines()
    TOOLS = [tool for tool in pitviper_res.keys() if tool != "DESeq2"]

    tool = widgets.Dropdown(options=TOOLS, value=TOOLS[0], description="Tool:")
    tool.observe(update_conditions, "value")

    conditions_list = list(pitviper_res[tool.value].keys())
    conditions = widgets.Dropdown(
        options=conditions_list, value=conditions_list[0], description="Condition:"
    )

    description = widgets.Text(
        value="My gene list", placeholder="Description", description="Description:"
    )
    bases = widgets.SelectMultiple(options=BASES)

    col_2 = widgets.ColorPicker(concise=False, description="Top color", value="blue")
    col_1 = widgets.ColorPicker(concise=False, description="Bottom color", value="red")
    plot_type = widgets.Dropdown(
        options=["Circle", "Bar"], value="Circle", description="Plot type:"
    )
    size = widgets.Dropdown(
        options=[5, 10, 20, 50, 100, 200, "max"], value=5, description="Size:"
    )
    fdr_cutoff = widgets.FloatSlider(
        min=0.0, max=1.0, step=0.01, value=0.05, description="FDR cut-off:"
    )
    score_cutoff = widgets.IntText(value=0, placeholder=0, description="Score cut-off:")
    button = widgets.Button(description="EnrichR!")

    display(
        widgets.VBox(
            [
                tool,
                conditions,
                description,
                bases,
                fdr_cutoff,
                score_cutoff,
                plot_type,
                size,
                col_2,
                col_1,
                button,
            ]
        )
    )

    def on_button_clicked(b):
        charts = []
        tool_res = pitviper_res[tool.value]
        treatment, baseline = conditions.value.split("_vs_")
        for base in bases.value:
            if tool.value == "MAGeCK_MLE":
                info = tool_res[conditions.value][
                    conditions.value + ".gene_summary.txt"
                ]
                info = info.loc[info[treatment + "|fdr"] < fdr_cutoff.value]
                if score_cutoff.value < 0:
                    info = info.loc[info[treatment + "|beta"] < score_cutoff.value]
                elif score_cutoff.value > 0:
                    info = info.loc[info[treatment + "|beta"] > score_cutoff.value]
                genes = info["Gene"].to_list()

            if tool.value == "MAGeCK_RRA":
                info = tool_res[conditions.value][
                    conditions.value + ".gene_summary.txt"
                ]
                info = info.loc[info["neg|fdr"] < fdr_cutoff.value]
                genes = info["id"]

            if tool.value == "BAGEL":
                info = tool_res[conditions.value][conditions.value + "_BAGEL_output.bf"]
                info = info.loc[info["BF"] > score_cutoff.value]
                genes = info["GENE"]

            if tool.value == "in_house_method":
                info = tool_res[conditions.value][
                    conditions.value + "_all-elements_in-house.txt"
                ]
                if score_cutoff.value > 0:
                    info = info.loc[info["score"] > score_cutoff.value]
                else:
                    info = info.loc[info["score"] < score_cutoff.value]
                genes = info["Gene"]

            if tool.value == "GSEA-like":
                info = tool_res[conditions.value][
                    conditions.value + "_all-elements_GSEA-like.txt"
                ]
                if score_cutoff.value > 0:
                    info = info.loc[info["NES"] > score_cutoff.value]
                elif score_cutoff.value < 0:
                    info = info.loc[info["NES"] < score_cutoff.value]
                info = info.loc[info["padj"] < fdr_cutoff.value]
                genes = info["pathway"]

            if tool.value == "CRISPhieRmix":
                info = tool_res[conditions.value][conditions.value + ".txt"]
                info = info.loc[info["locfdr"] < fdr_cutoff.value]
                if score_cutoff.value > 0:
                    info = info.loc[
                        info["top_3_mean_log2FoldChange"] > score_cutoff.value
                    ]
                elif score_cutoff.value <= 0:
                    info = info.loc[
                        info["top_3_mean_log2FoldChange"] < score_cutoff.value
                    ]
                genes = info["gene"]

            enrichr_res = getEnrichrResults(genes, description.value, base)
            table = createEnrichrTable(enrichr_res)
            if plot_type.value == "Bar":
                chart = enrichmentBarPlot(
                    table, size.value, description.value, col_1.value, col_2.value, base
                )
            else:
                chart = enrichmentCirclePlot(
                    table, size.value, description.value, col_1.value, col_2.value, base
                )
            charts.append(chart)
        for chart in charts:
            display(chart)

    button.on_click(on_button_clicked)


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
        # return "This module is available only if genes symbol are used."

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
            info = tool_res[conditions.value][conditions.value + "_BAGEL_output.bf"]
            info = info.loc[info["BF"] > float(score_cutoff.value)]
            genes = info["GENE"]

        if tool.value == "in_house_method":
            info = tool_res[conditions.value][
                conditions.value + "_all-elements_in-house.txt"
            ]
            if score_cutoff.value > 0:
                info = info.loc[info["score"] > score_cutoff.value]
            else:
                info = info.loc[info["score"] < score_cutoff.value]
            genes = info["Gene"]

        if tool.value == "GSEA-like":
            info = tool_res[conditions.value][
                conditions.value + "_all-elements_GSEA-like.txt"
            ]
            if float(score_cutoff.value) > 0:
                info = info.loc[info["NES"] > float(score_cutoff.value)]
            elif float(score_cutoff.value) < 0:
                info = info.loc[info["NES"] < float(score_cutoff.value)]
            info = info.loc[info["padj"] < fdr_cutoff.value]
            genes = info["pathway"]

        print("Size (gene set):", len(genes))
        link = "http://genemania.org/search/homo-sapiens/" + "/".join(genes)
        print(link)

    BASES = open("workflow/notebooks/enrichr_list.txt", "r").readlines()
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


def ranking(treatment, control, token, tools_available, params):
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
        score = params["BAGEL"]["score"]
        greater = params["BAGEL"]["greater"]
        bagel = tools_available["BAGEL"][comparison][comparison + "_BAGEL_output.bf"]
        if greater:
            bagel = bagel[(bagel["BF"] > score)]
        else:
            bagel = bagel[(bagel["BF"] < score)]
        bagel["default_rank"] = bagel["BF"].rank(method="dense", ascending=False).copy()
        bagel = bagel[["GENE", "default_rank"]].rename(
            columns={"GENE": "id", "default_rank": "bagel_rank"}
        )
        bagel_genes = list(bagel.id)
        tool_results["BAGEL"] = bagel_genes
        tool_genes.append(bagel_genes)

    if params["in_house_method"]["on"]:
        in_house = tools_available["in_house_method"][comparison][
            comparison + "_all-elements_in-house.txt"
        ]
        if params["in_house_method"]["direction"] == "Negative":
            in_house = in_house.loc[in_house.category == "down"]
        elif params["in_house_method"]["direction"] == "Positive":
            in_house = in_house.loc[in_house.category == "up"]
        in_house["default_rank"] = in_house["score"].rank(method="first").copy()
        in_house = in_house[["Gene", "default_rank"]].rename(
            columns={"Gene": "id", "default_rank": "in_house_rank"}
        )
        in_house_genes = list(in_house.id)
        tool_results["In-house"] = in_house_genes
        tool_genes.append(in_house_genes)

    if params["GSEA_like"]["on"]:
        score = params["GSEA_like"]["score"]
        fdr = params["GSEA_like"]["fdr"]
        greater = params["GSEA_like"]["greater"]
        gsea = tools_available["GSEA-like"][comparison][
            comparison + "_all-elements_GSEA-like.txt"
        ]
        if not greater:
            gsea = gsea[(gsea["NES"] < score) & (gsea["padj"] < fdr)]
        else:
            gsea = gsea[(gsea["NES"] > score) & (gsea["padj"] < fdr)]
        gsea["default_rank"] = gsea["NES"].rank(method="dense").copy()
        gsea = gsea[["pathway", "default_rank"]].rename(
            columns={"pathway": "id", "default_rank": "gsea_rank"}
        )
        gsea_genes = list(gsea.id)
        tool_results["GSEA-like"] = gsea_genes
        tool_genes.append(gsea_genes)

    if params["CRISPhieRmix"]["on"]:
        score = params["CRISPhieRmix"]["top_3_mean_log2FoldChange"]
        fdr = params["CRISPhieRmix"]["fdr"]
        greater = params["CRISPhieRmix"]["greater"]
        crisphie = tools_available["CRISPhieRmix"][comparison][comparison + ".txt"]
        if not greater:
            crisphie = crisphie[
                (crisphie["top_3_mean_log2FoldChange"] < score)
                & (crisphie["locfdr"] < fdr)
            ]
        else:
            crisphie = crisphie[
                (crisphie["top_3_mean_log2FoldChange"] > score)
                & (crisphie["locfdr"] < fdr)
            ]
        crisphie["default_rank"] = crisphie["FDR"].rank(method="dense").copy()
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
        bagel = tools_available["BAGEL"][comparison][comparison + "_BAGEL_output.bf"]
        bagel["default_rank"] = bagel["BF"].rank(method="dense", ascending=False).copy()
        bagel = bagel[["GENE", "default_rank"]].rename(
            columns={"GENE": "id", "default_rank": "bagel_rank"}
        )
        pdList.append(bagel)

    if params["in_house_method"]["on"]:
        in_house = tools_available["in_house_method"][comparison][
            comparison + "_all-elements_in-house.txt"
        ]
        in_house["default_rank"] = in_house["score"].rank(method="first").copy()
        in_house = in_house[["Gene", "default_rank"]].rename(
            columns={"Gene": "id", "default_rank": "in_house_rank"}
        )
        pdList.append(in_house)

    if params["GSEA_like"]["on"]:
        gsea = tools_available["GSEA-like"][comparison][
            comparison + "_all-elements_GSEA-like.txt"
        ]
        gsea["default_rank"] = gsea["NES"].rank(method="dense").copy()
        gsea = gsea[["pathway", "default_rank"]].rename(
            columns={"pathway": "id", "default_rank": "gsea_rank"}
        )
        pdList.append(gsea)

    if params["CRISPhieRmix"]["on"]:
        crisphie = tools_available["CRISPhieRmix"][comparison][comparison + ".txt"]
        crisphie["default_rank"] = crisphie["FDR"].rank(method="dense").copy()
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
    venn_lib = importr("venn")
    grdevices = importr("grDevices")
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_occurences = ro.conversion.py2rpy(occurences)

    with rpy2.robjects.lib.grdevices.render_to_bytesio(
        grdevices.png, width=1024, height=896, res=150
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
            "top_3_mean_log2FoldChange": 0,
            "greater": False,
        },
        "in_house_method": {"on": False, "direction": "Negative"},
        "GSEA_like": {"on": False, "fdr": 0.25, "score": 0, "greater": False},
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
        params["CRISPhieRmix"]["top_3_mean_log2FoldChange"] = change["new"]

    ### In-house
    def in_house_method_direction_update(change):
        params["in_house_method"]["direction"] = change["new"]

    ### GSEA-like
    def GSEA_like_order_update(change):
        if change["new"] == "Greater than score":
            params["GSEA_like"]["greater"] = True
        else:
            params["GSEA_like"]["greater"] = False

    def GSEA_like_fdr_update(change):
        params["GSEA_like"]["fdr"] = change["new"]

    def GSEA_like_score_update(change):
        params["GSEA_like"]["score"] = change["new"]

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
        bagel_score = widgets.FloatText(value=0, description="Score cut-off:")
        bagel_text = widgets.HTML(value="<b>BAGEL</b>:")
        bagel_order = widgets.ToggleButtons(
            options=["Greater than score", "Lower than score"], description="Selection:"
        )
        display(bagel_text)
        bagel_box = widgets.HBox([bagel_score, bagel_order])
        display(bagel_box)
        bagel_order.observe(bagel_order_update, "value")
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
    if "in_house_method" in tools_selected:
        params["in_house_method"]["on"] = True
        in_house_method_direction = widgets.ToggleButtons(
            options=["Negative", "Positive"], description="Direction:"
        )
        in_house_method_text = widgets.HTML(value="<b>In-house method</b>:")
        display(in_house_method_text)
        in_house_method_box = widgets.HBox([in_house_method_direction])
        display(in_house_method_box)
        in_house_method_direction.observe(in_house_method_direction_update, "value")
    if "GSEA-like" in tools_selected:
        params["GSEA_like"]["on"] = True
        GSEA_like_fdr = widgets.FloatSlider(
            min=0.0, max=1.0, step=0.01, value=0.25, description="FDR:"
        )
        GSEA_like_score = widgets.FloatText(value=0, description="Score cut-off:")
        GSEA_like_text = widgets.HTML(value="<b>GSEA-like</b>:")
        GSEA_like_order = widgets.ToggleButtons(
            options=["Lower than score", "Greater than score"], description="Selection:"
        )
        display(GSEA_like_text)
        GSEA_like_box = widgets.HBox([GSEA_like_fdr, GSEA_like_score, GSEA_like_order])
        display(GSEA_like_box)
        GSEA_like_order.observe(GSEA_like_order_update, "value")
        GSEA_like_score.observe(GSEA_like_score_update, "value")
        GSEA_like_fdr.observe(GSEA_like_fdr_update, "value")


def disable_widgets(token):
    config = "./config/%s.yaml" % token
    content = open_yaml(config)
    disabled = False
    if content["screen_type"] == "not_gene":
        disabled = True
    return disabled


def multiple_tools_results(tools_available, token):
    TOOLS = [tool for tool in tools_available.keys() if not tool in ["DESeq2"]]

    # Define widgets's options
    conditions_options = tools_available[TOOLS[0]].keys()
    tools_options = TOOLS

    # Define widgets
    conditions_widget = widgets.Dropdown(
        options=conditions_options, description="Conditions:"
    )
    tools_widget = widgets.SelectMultiple(options=tools_options, description="Tool:")
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
        fig = ff.create_dendrogram(data_array, orientation="bottom")  # , labels=labels)
        for i in range(len(fig["data"])):
            fig["data"][i]["yaxis"] = "y2"

        # Create Side Dendrogram
        dendro_side = ff.create_dendrogram(
            data_array.transpose(), orientation="right"
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

        # # Edit yaxis
        fig.update_layout(
            font={"size": 10},
            yaxis={
                "domain": [0, 0.85],
                "ticktext": cols[indicies_yaxis],
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
        description="Venn diagram", style=dict(button_color="#707070")
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
            if params["BAGEL"]["greater"]:
                word = "greater"
            else:
                word = "less"
            display(
                HTML(
                    """<p style="color:black;padding: 0.5em;"><b>BAGEL</b>: Bayesian factor threshold = %s, keep elements with Bayesian factor %s than threshold.</p>"""
                    % (params["BAGEL"]["score"], word)
                )
            )
        if "CRISPhieRmix" in tools_widget.value:
            if params["CRISPhieRmix"]["greater"]:
                word = "greater"
            else:
                word = "less"
            display(
                HTML(
                    """<p style="color:black;padding: 0.5em;"><b>CRISPhieRmix</b>: FDR < %s, Top-3 LogFoldChange Average threshold = %s, keep elements with Top-3 LogFoldChange Average %s than threshold.</p>"""
                    % (
                        params["CRISPhieRmix"]["fdr"],
                        params["CRISPhieRmix"]["top_3_mean_log2FoldChange"],
                        word,
                    )
                )
            )
        if "GSEA-like" in tools_widget.value:
            if params["GSEA_like"]["greater"]:
                word = "greater"
            else:
                word = "less"
            display(
                HTML(
                    """<p style="color:black;padding: 0.5em;"><b>GSEA-like (fGSEA)</b>: FDR < %s, NES threshold = %s, keep elements with NES %s than threshold.</p>"""
                    % (params["GSEA_like"]["fdr"], params["GSEA_like"]["score"], word)
                )
            )
        if "in_house_method" in tools_widget.value:
            display(
                HTML(
                    """<p style="color:black;padding: 0.5em;"><b>Filtering method</b>: %s selection.</p>"""
                    % (params["in_house_method"]["direction"])
                )
            )

    def venn_button_clicked(b):
        treatment, control = conditions_widget.value.split("_vs_")
        ranks, occurences = ranking(treatment, control, token, tools_available, params)
        if selection_widgets.value == "Intersection":
            df = pd.DataFrame(
                occurences.eq(occurences.iloc[:, 0], axis=0).all(1),
                columns=["intersection"],
            )
            genes_list = df.loc[df.intersection == True].index
        else:
            df = pd.DataFrame(
                occurences.eq(occurences.iloc[:, 0], axis=0).any(1), columns=["union"]
            )
            genes_list = df.loc[df.union == True].index
        display(
            HTML(
                """<p style="color:white;font-weight: bold;background-color: orange;padding: 0.5em;">Venn diagram: %s</p>"""
                % selection_widgets.value
            )
        )
        show_parameters(params)
        plot_venn(occurences)
        print("Genes at %s of all methods:" % selection_widgets.value)
        for gene in genes_list:
            print(gene)

    def rra_button_clicked(b):
        treatment, control = conditions_widget.value.split("_vs_")
        ranks, occurences = ranking(treatment, control, token, tools_available, params)
        display(
            HTML(
                """<p style="color:white;font-weight: bold;background-color: green;padding: 0.5em;">RRA results</p>"""
            )
        )
        show_parameters(params)
        run_rra(ranks)

    def genemania_button_clicked(b):
        treatment, control = conditions_widget.value.split("_vs_")
        ranks, occurences = ranking(treatment, control, token, tools_available, params)
        if selection_widgets.value == "Intersection":
            df = pd.DataFrame(
                occurences.eq(occurences.iloc[:, 0], axis=0).all(1),
                columns=["intersection"],
            )
            genes_list = df.loc[df.intersection == True].index
        else:
            df = pd.DataFrame(
                occurences.eq(occurences.iloc[:, 0], axis=0).any(1), columns=["union"]
            )
            genes_list = df.loc[df.union == True].index
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
            b, genes, bases, size, plot_type, col_2, col_1, description
        ):
            charts = []
            title = description.value + " (%s)" % selection_widgets.value
            for base in bases.value:
                enrichr_res = getEnrichrResults(genes, description.value, base)
                table = createEnrichrTable(enrichr_res)
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

        display(
            HTML(
                """<p style="color:white;font-weight: bold;background-color: purple;padding: 0.5em;">EnrichR for genes at %s</p>"""
                % selection_widgets.value
            )
        )
        show_parameters(params)
        treatment, control = conditions_widget.value.split("_vs_")
        ranks, occurences = ranking(treatment, control, token, tools_available, params)
        if selection_widgets.value == "Intersection":
            df = pd.DataFrame(
                occurences.eq(occurences.iloc[:, 0], axis=0).all(1),
                columns=["intersection"],
            )
            genes_list = df.loc[df.intersection == True].index
        else:
            df = pd.DataFrame(
                occurences.eq(occurences.iloc[:, 0], axis=0).any(1), columns=["union"]
            )
            genes_list = df.loc[df.union == True].index
        BASES = open("workflow/notebooks/enrichr_list.txt", "r").readlines()
        bases = widgets.SelectMultiple(options=BASES, description="Gene sets:", rows=10)
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
            value="My gene list", placeholder="Description", description="Description:"
        )
        button_enrichr = widgets.Button(description="EnrichR!")

        display(
            widgets.VBox(
                [description, bases, plot_type, size, col_2, col_1, button_enrichr]
            )
        )
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
            )
        )

    def depmap_button_clicked(b):
        display(
            HTML(
                """<p style="color:white;font-weight: bold;background-color: #A52A2A;padding: 0.5em;">depmap vizualisation module</p>"""
            )
        )
        data_types_widget = widgets.RadioButtons(
            options=["crispr", "proteomic", "rnai", "tpm", "mutations"],
            value="crispr",
            description="Data type:",
            disabled=False,
        )

        def download_depmap_file(data_type, release):
            target_file = "resources/depmap/%s_%s.txt" % (release, data_type)
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

        def getRelease():
            depmap_release = depmap.depmap_release()
            return str(depmap_release).rstrip()[5:-1]

        def depmap_query_button_clicked(b):
            depmap_release = getRelease()
            save_path = "resources/depmap/%s_%s.txt" % (
                depmap_release,
                data_types_widget.value,
            )
            Path("resources/depmap/").mkdir(parents=True, exist_ok=True)
            treatment, control = conditions_widget.value.split("_vs_")
            ranks, occurences = ranking(
                treatment, control, token, tools_available, params
            )
            if download_depmap_file(data_types_widget.value, depmap_release):
                print("This step can take some time.")
                print("Querying: %s..." % data_types_widget.value)
                eh = experimentHub.ExperimentHub()
                base_package = importr("base")
                dplyr = importr("dplyr")
                tidyr = importr("tidyr")
                readr = importr("readr")
                if data_types_widget.value == "rnai":
                    depmap_data = tidyr.drop_na(
                        dplyr.select(
                            depmap.depmap_rnai(), "depmap_id", "gene_name", "dependency"
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
                    "resources/depmap/%s_metadata.txt" % depmap_release
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
                        "resources/depmap/%s_metadata.txt" % depmap_release,
                        row_names=False,
                        quote=False,
                        sep="\t",
                    )
                else:
                    print("Import metadata...")
                    depmap_metadata = readr.read_delim(
                        "resources/depmap/%s_metadata.txt" % depmap_release, delim="\t"
                    )
                depmap_data = base_package.merge(
                    depmap_data, depmap_metadata, by="depmap_id"
                )
                print("Saving %s" % save_path)
                utils.write_table(
                    depmap_data, save_path, row_names=False, quote=False, sep="\t"
                )
            print("Opening %s" % save_path)
            data = pd.read_table(save_path, sep="\t")

            tissues_init = list(set(data.cell_line))
            tissues = [
                "_".join(str(tissu).split("_")[1:])
                for tissu in tissues_init
                if not str(tissu) in ["nan", ""]
            ]
            tissues = list(set(tissues))
            tissues.insert(0, "All")

            primary_diseases = list(set(data.primary_disease))
            primary_diseases.insert(0, "All")

            cell_lines_init = list(set(data.cell_line_name))
            cell_lines = [
                tissu for tissu in cell_lines_init if not str(tissu) in ["nan", ""]
            ]
            cell_lines = natural_sort(cell_lines)
            cell_lines.insert(0, "All")

            tissues_widget = widgets.SelectMultiple(
                options=tissues, value=["All"], description="Tissu:"
            )
            primary_diseases_widget = widgets.SelectMultiple(
                options=primary_diseases, value=["All"], description="Primary tissu:"
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
                    essential_genes = table[table.gene_name.isin(genes_list)][columns]
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
                        alt.Chart(essential_genes, title="depmap deleterious mutations")
                        .mark_rect()
                        .encode(
                            x=alt.X("cell_line_name", axis=alt.Axis(title="Cell line")),
                            y=alt.Y("gene_name", axis=alt.Axis(title="Gene name")),
                            color=alt.Color(
                                "primary_disease", scale=alt.Scale(scheme="tableau20")
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
                    ]  # .pivot_table(index='gene_name', columns='cell_line', values=variable).dropna()
                    plot_interactive_heatmap(essential_genes, variable)

            button = widgets.Button(description="Run!")
            button.on_click(tissu_selection_button_clicked)
            display(tissues_widget, primary_diseases_widget, cell_lines_widget, button)

        depmap_query_button = widgets.Button(description="Query!")
        depmap_query_button.on_click(depmap_query_button_clicked)
        display(data_types_widget, depmap_query_button)

    venn_button.on_click(venn_button_clicked)
    rra_button.on_click(rra_button_clicked)
    genemania_button.on_click(genemania_button_clicked)
    enrichr_button.on_click(enrichr_button_clicked)
    depmap_button.on_click(depmap_button_clicked)

    def ranking_button_clicked(event):
        treatment, control = conditions_widget.value.split("_vs_")
        ranks, occurences = ranking(treatment, control, token, tools_available, params)
        ranks_path = Path("results/%s/Integration/ranks.csv" % token)
        ranks_path.parent.mkdir(parents=True, exist_ok=True)
        ranks.to_csv(ranks_path)
        print("Ranks table writed in %s" % ranks_path)
        occurences_path = Path("results/%s/Integration/occurences.csv" % token)
        occurences_path.parent.mkdir(parents=True, exist_ok=True)
        occurences.to_csv(occurences_path)
        print("Occurences table writed in %s" % occurences_path)

    ranking_button.on_click(ranking_button_clicked)


def plot_chart(
    data,
    tool,
    condition,
    file,
    x,
    y,
    method,
    column_filter="",
    cutoff=0,
    greater=False,
    element="",
    column_element="",
    elem_color="",
    pass_color="",
    fail_color="",
    category_column="",
    color_scheme="",
    reverse="",
):
    source = data[tool][condition][file].copy()
    show_plot = True
    if method == "Cut-off":
        if source[column_filter].dtype in [np.float64, np.int64]:
            elements = element.split(",")
            if not greater:
                pass_value = "%s < %s" % (column_filter, cutoff)
                fail_value = "%s >= %s" % (column_filter, cutoff)
                source.loc[source[column_filter] < cutoff, "filter"] = pass_value
                source.loc[source[column_filter] >= cutoff, "filter"] = fail_value
            else:
                pass_value = "%s >= %s" % (column_filter, cutoff)
                fail_value = "%s < %s" % (column_filter, cutoff)
                source.loc[source[column_filter] >= cutoff, "filter"] = pass_value
                source.loc[source[column_filter] < cutoff, "filter"] = fail_value
            source["filter"] = np.where(
                source[column_element].isin(elements), element, source["filter"]
            )
            domain = [pass_value, fail_value, element]
            range_ = [pass_color, fail_color, elem_color]
        else:
            source["filter"] = "pass"
            show_plot = False
            print("Choose a column of integer or float type.")
    if show_plot:
        if method == "Cut-off":
            chart = (
                alt.Chart(source)
                .transform_calculate(
                    order="{'%s': 2, '%s': 1, '%s': 0}[datum.filter]"
                    % (element, pass_value, fail_value)
                )
                .mark_circle()
                .encode(
                    x=x,
                    y=y,
                    tooltip=list(source.columns),
                    color=alt.Color(
                        "filter",
                        scale=alt.Scale(domain=domain, range=range_),
                        legend=alt.Legend(title="Pass filter:"),
                    ),
                    order="order:Q",
                )
                .interactive()
            )
            display(chart)
        else:
            chart = (
                alt.Chart(source)
                .mark_circle()
                .encode(
                    x=x,
                    y=y,
                    tooltip=list(source.columns),
                    color=alt.Color(
                        category_column,
                        legend=alt.Legend(title="%s:" % category_column),
                        scale=alt.Scale(scheme=color_scheme, reverse=reverse),
                    ),
                )
                .interactive()
            )
            print(chart)


def call_form(tools_available):
    SCHEMES = ["reds", "accent", "redblue", "rainbow"]
    tools = list(tools_available.keys())
    display(widgets.HTML(value="First, choose a <b>tool</b> to browse results:"))

    @interact(
        tool=widgets.Dropdown(
            options=tools, value=tools[0], description="Tool:", disabled=False
        )
    )
    def form(tool):
        display(
            widgets.HTML(value="Choose a <b>condition</b> available for this tool:")
        )
        conditions = list(tools_available[tool].keys())

        @interact(
            condition=widgets.Dropdown(
                options=conditions,
                value=conditions[0],
                description="Condition:",
                disabled=False,
            )
        )
        def form(condition):
            display(widgets.HTML(value="Choose a <b>file</b>:"))
            files = list(tools_available[tool][condition].keys())

            @interact(
                file=widgets.Dropdown(
                    options=files, value=files[0], description="File:", disabled=False
                )
            )
            def form(file):
                display(
                    widgets.HTML(
                        value="Choose which columns to use as <b>x</b> and <b>y</b> axis:"
                    )
                )
                columns = tools_available[tool][condition][file].columns

                @interact(
                    x=widgets.Dropdown(
                        options=columns,
                        value=columns[0],
                        description="x:",
                        disabled=False,
                    ),
                    y=widgets.Dropdown(
                        options=columns,
                        value=columns[1],
                        description="y:",
                        disabled=False,
                    ),
                )
                def form(x, y):
                    display(
                        widgets.HTML(
                            value="Results can be colored using two <b>methods</b>.</br> <ul><li>Cut-off: Choose a <b>column</b> with numerical values on which you want to apply a <b>cut-off</b></li><li>Category: choose a column and a color scheme.</li></ul>"
                        )
                    )

                    @interact(
                        method=widgets.ToggleButtons(
                            options=["Cut-off", "Category"],
                            value="Cut-off",
                            description="Method:",
                        )
                    )
                    def form(method):
                        if method == "Cut-off":
                            display(
                                widgets.HTML(
                                    value="Now choose a column with numerical data and a set a <b>cut-off</b> value."
                                )
                            )

                            @interact(
                                column_filter=widgets.Dropdown(
                                    options=columns,
                                    value=columns[0],
                                    description="Column:",
                                    disabled=False,
                                ),
                                cutoff=widgets.FloatText(
                                    value=0.05, description="Cut-off:"
                                ),
                                greater=widgets.Checkbox(
                                    value=False,
                                    description="Greater?",
                                    disabled=False,
                                    indent=True,
                                ),
                                pass_color=widgets.ColorPicker(
                                    concise=False,
                                    description="Pass color:",
                                    value="red",
                                    disabled=False,
                                ),
                                fail_color=widgets.ColorPicker(
                                    concise=False,
                                    description="Fail color:",
                                    value="gray",
                                    disabled=False,
                                ),
                                column_element=widgets.Dropdown(
                                    options=columns,
                                    value=columns[0],
                                    description="Column:",
                                    disabled=False,
                                ),
                            )
                            def form(
                                column_filter,
                                cutoff,
                                greater,
                                pass_color,
                                fail_color,
                                column_element,
                            ):
                                display(
                                    widgets.HTML(
                                        value="It is possible to highlight an <b>element</b> based on a column and value from that column."
                                    )
                                )
                                first_element = str(
                                    tools_available[tool][condition][file][
                                        column_element
                                    ][0]
                                )

                                @interact(
                                    element=widgets.Text(
                                        value=first_element,
                                        placeholder="Element:",
                                        description="Element:",
                                        disabled=False,
                                    ),
                                    elem_color=widgets.ColorPicker(
                                        concise=False,
                                        description="Element color:",
                                        value="blue",
                                        disabled=False,
                                    ),
                                )
                                def form(element, elem_color):
                                    plot_chart(
                                        tools_available,
                                        tool,
                                        condition,
                                        file,
                                        x,
                                        y,
                                        method,
                                        column_filter,
                                        cutoff,
                                        greater,
                                        element,
                                        column_element,
                                        elem_color,
                                        pass_color,
                                        fail_color,
                                    )

                        elif method == "Category":

                            @interact(
                                category_column=widgets.Dropdown(
                                    options=columns,
                                    value=columns[0],
                                    description="Category:",
                                    disabled=False,
                                ),
                                reverse=widgets.Checkbox(
                                    value=False,
                                    description="Reverse?",
                                    disabled=False,
                                    indent=True,
                                ),
                                color_scheme=widgets.Dropdown(
                                    options=SCHEMES,
                                    value=SCHEMES[0],
                                    description="Scheme:",
                                    disabled=False,
                                ),
                            )
                            def form(category_column, reverse, color_scheme):
                                plot_chart(
                                    tools_available,
                                    tool,
                                    condition,
                                    file,
                                    x,
                                    y,
                                    method,
                                    category_column=category_column,
                                    color_scheme=color_scheme,
                                    reverse=reverse,
                                )
