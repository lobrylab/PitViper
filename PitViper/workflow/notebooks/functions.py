import os
import pandas as pd
pd.options.mode.chained_assignment = None
import re
import altair as alt
from ipywidgets import interact, interactive, fixed, interact_manual, Button, HBox, VBox
import ipywidgets as widgets
import requests
import json
from pathlib import Path
import natsort as ns
import numpy as np
import yaml
import os.path
from os import path
from os import listdir
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import decomposition
from sklearn import datasets
from functools import reduce
from IPython.display import display, clear_output
from clustergrammer2 import net, Network, CGM2
import rpy2
import IPython
import math
import functools
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as robjects
from rpy2.robjects.conversion import rpy2py

import gc
######################################################
dplyr = importr("dplyr")
tibble = importr("tibble")
stringr = importr("stringr")
depmap = importr("depmap")
experimentHub = importr("ExperimentHub")
utils = importr('utils')

from rpy2.robjects.lib.dplyr import DataFrame
import rpy2.ipython.html
rpy2.ipython.html.init_printing()
######################################################

alt.renderers.enable('html')


def working_directory_update(output):
    switch = True
    while os.path.basename(os.getcwd()) != "PitViper":
        if switch:
            switch = False
            os.system("cd ../../")
        else:
            os.system("cd ../")

    print('Working directory: ', os.getcwd())

    with open(output, "w") as out:
        print("Notebook was runned.", file=out)


def open_yaml(yml):
    """Open a YAML file and return it's content."""
    with open(yml, "r") as stream:
        try:
            content = yaml.safe_load(stream)
            return content
        except yaml.YAMLError as exc:
            print(exc)


def import_results(token):
    print('Token: %s\n' % token)
    config = "./config/%s.yaml" % token
    print('Config file used: %s' % config)
    content = open_yaml(config)

    tools = ["DESeq2"]

    if content["mageck_mle_activate"] == 'True':
        tools.append("MAGeCK_MLE")
    if content["mageck_rra_activate"] == 'True':
        tools.append("MAGeCK_RRA")
    if content["bagel_activate"] == 'True':
        tools.append("BAGEL")
    if content["crisphiermix_activate"] == 'True':
        tools.append("CRISPhieRmix")
    if content["filtering_activate"] == 'True':
        tools.append("in_house_method")
    if content["gsea_activate"] == 'True':
        tools.append("GSEA-like")

    results_directory = "results/%s/" % token
    print("Results directory: %s \n" % results_directory)

    tools_available = {}
    print('Tools available:')
    for tool in tools:
        if tool in os.listdir(results_directory):
            print("\t-%s" % tool)
            tools_available[tool] = {}

    for tool in tools_available:
        print("- Process %s results..." % tool)
        for comparison in os.listdir(results_directory + tool):
            tools_available[tool][comparison] = {}
            for file in os.listdir(os.path.join(results_directory, tool, comparison)):
                if file.endswith('.txt') or file.endswith('.bf'):
                    if tool in ["CRISPhieRmix"]:
                        sep = ","
                    else:
                        sep = "\t"
                    tools_available[tool][comparison][file] = pd.read_csv(os.path.join(results_directory, tool, comparison, file), sep = sep)
    add_columns(tools_available)
    return (results_directory, tools_available)


def add_columns(tools_available):
    if 'MAGeCK_MLE' in tools_available:
        for condition in tools_available['MAGeCK_MLE']:
            treatment = condition.split('_vs_')[0]
            for file_suffixe in ['%s.gene_summary.txt', '%s.sgrna_summary.txt']:
                table = tools_available['MAGeCK_MLE'][condition][file_suffixe % condition]
                if (not 'log10(invFDR)' in list(table.columns))  and ('%s|fdr' % treatment in list(table.columns)):
                    array = table['%s|fdr' % treatment].values
                    min_fdr = np.min(array[np.nonzero(array)])
                    table['%s|fdr_nozero' % treatment] = table['%s|fdr' % treatment].replace(0, min_fdr)
                    min_fdr = table['%s|fdr' % treatment].values
                    table['log10(invFDR)'] = - np.log2(table['%s|fdr_nozero' % treatment])
    if 'MAGeCK_RRA' in tools_available:
        for condition in tools_available['MAGeCK_RRA']:
            treatment = condition.split('_vs_')[0]
            for file_suffixe in ['%s.gene_summary.txt', '%s.sgrna_summary.txt']:
                table = tools_available['MAGeCK_RRA'][condition][file_suffixe % condition]
                for direction in ['neg', 'pos']:
                    if not '%s|log10(invFDR)' % direction in list(table.columns) and ('%s|fdr' % direction in list(table.columns)):
                        array = table['%s|fdr' % direction].values
                        min_fdr = np.min(array[np.nonzero(array)])
                        table['%s|fdr_nozero' % direction] = table['%s|fdr' % direction].replace(0, min_fdr)
                        min_fdr = table['%s|fdr' % direction].values
                        table['%s|log10(invFDR)' % direction] = - np.log2(table['%s|fdr_nozero' % direction])
    if 'GSEA-like' in tools_available:
        for condition in tools_available['GSEA-like']:
            treatment = condition.split('_vs_')[0]
            for file_suffixe in ['%s_all-elements_GSEA-like.txt']:
                table = tools_available['GSEA-like'][condition][file_suffixe % condition]
                if (not 'log10(invPadj)' in list(table.columns))  and ('padj' in list(table.columns)):
                    array = table['padj'].values
                    min_fdr = np.min(array[np.nonzero(array)])
                    table['padj_nozero'] = table['padj'].replace(0, min_fdr)
                    min_fdr = table['padj'].values
                    table['log10(invPadj)'] = - np.log2(table['padj'])


def show_mapping_qc(token):
    path_qc = "./resources/%s/screen.countsummary.txt" % token
    if not path.exists(path_qc):
        print("No mapping QC file to show.")
        return 0
    table = pd.read_csv(path_qc, sep='\t')
    table = table[['Label', 'Reads', 'Mapped', 'Percentage', 'Zerocounts', 'GiniIndex']]
    table['Label'] = pd.Categorical(table['Label'], ordered=True, categories= ns.natsorted(table['Label'].unique()))
    table = table.sort_values('Label')

    def color_low_mapping_red(val):
        """
        Takes a scalar and returns a string with
        the css property `'color: red'` for negative
        strings, black otherwise.
        """
        color = 'red' if float(val) < 0.6 else 'green'
        return 'color: %s' % color

    def color_high_gini_red(val):
        """
        Takes a scalar and returns a string with
        the css property `'color: red'` for negative
        strings, black otherwise.
        """
        color = 'red' if val > 0.35 else 'green'
        return 'color: %s' % color

    s = table.style.\
        applymap(color_low_mapping_red, subset=['Percentage']).\
        applymap(color_high_gini_red, subset=['GiniIndex'])
    return s


def show_read_count_distribution(token, width=800, height=400):
    config = "./config/%s.yaml" % token
    content = open_yaml(config)
    path_qc = content['normalized_count_table']
    if not path.exists(path_qc):
        print("No count file to show.")
        return 0
    table = pd.read_csv(path_qc, sep='\t')

    table.iloc[:, 2:] = table.iloc[:, 2:] +1
    table.iloc[:, 2:] = table.iloc[:, 2:].apply(np.log2)

    chart = alt.Chart(table).transform_fold(
        list(table.columns[2:]),
        as_ = ['Measurement_type', 'counts']
    ).transform_density(
        density='counts',
        bandwidth=0.3,
        groupby=['Measurement_type'],
        extent= [0, 20],
        counts = True,
        steps=200
    ).mark_line().encode(
        alt.X('value:Q', axis=alt.Axis(title='log2(read count)')),
        alt.Y('density:Q'),
        alt.Color('Measurement_type:N'),
        tooltip = ['Measurement_type:N', 'value:Q', 'density:Q'],
    ).properties(width=width, height=height)

    return chart



def natural_sort(l):
    """Function for natural sorting of list."""
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def CRISPhieRmix_data(comparison = "", control = "", tool = "", results_directory = "", tools_available = ""):
    tables_list = []
    for _comparison in os.listdir(os.path.join(results_directory, tool)):
        if _comparison.split("_vs_")[-1] == control:
            keys_list = list(tools_available[tool][_comparison].keys())
            data = tools_available[tool][_comparison][keys_list[0]]
            trt = _comparison.split("_vs_")[0]
            data['condition'] = trt
            tables_list.append(data)
        if _comparison == comparison:
            data = tools_available[tool][_comparison]["%s.txt" % comparison]
            tables_list.append(data)
    result = pd.concat(tables_list)
    return result



def CRISPhieRmix_results(results_directory, tools_available):
    tool = "CRISPhieRmix"
    if not tool in tools_available.keys():
        return "CRISPhieRmix wasn't activated."
    comparisons_list = os.listdir(os.path.join(results_directory, tool))
    ctrs = list(map(lambda x: x.split("_vs_")[1], list(tools_available[tool].keys())))
    @interact(control=widgets.Dropdown(options=set(ctrs), description='Control:', disabled=False),
              fdr_cutoff=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05),
              gene=widgets.Text(value='', placeholder='Feature to show...', description='Feature:'))
    def CRISPhieRmix_results(control, fdr_cutoff, gene):
        tool = "CRISPhieRmix"
        result = CRISPhieRmix_data(comparison = "", control = control, tool = tool, results_directory=results_directory, tools_available=tools_available)
        result.loc[result['locfdr'] < fdr_cutoff, 'significant'] = 'Yes'
        result.loc[result['locfdr'] >= fdr_cutoff, 'significant'] = 'No'
        new_row = {'gene':gene, 'condition':control, 'significant': 'Baseline', 'locfdr':1, 'score': 0}
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.gene == gene]
        domain = ['Yes', 'No', 'Baseline']
        range_ = ['red', 'grey', 'black']



        def on_button_clicked(b):
            sort_cols = natural_sort(list(res.condition.values))
            plot = alt.Chart(res).mark_circle(size=60).mark_point(
                filled=True,
                size=100,
                ).encode(
                        y='locfdr',
                        x=alt.X('condition:O', sort=sort_cols),
                        color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
                        tooltip=['gene', 'locfdr', 'significant', 'score'],
                ).properties(
                        title=gene + " locfdr versus baseline (CRISPhieRmix)",
                        width=100
                )
            with output:
                display(plot)
                return(plot)

        button = widgets.Button(description="Show plot")
        output = widgets.Output()
        display(button, output)
        button.on_click(on_button_clicked)



def GSEA_like_data(comparison = "", control = "", tool = "", results_directory = "", tools_available = ""):
    tables_list = []
    for _comparison in os.listdir(os.path.join(results_directory, tool)):
        if _comparison.split("_vs_")[-1] == control:
            keys_list = list(tools_available[tool][_comparison].keys())
            for key in keys_list:
                if key.endswith("_all-elements_GSEA-like.txt"):
                    break
            data = tools_available[tool][_comparison][key]
            trt = _comparison.split("_vs_")[0]
            data['condition'] = trt
            tables_list.append(data)
        if _comparison == comparison:
            data = tools_available[tool][_comparison]["%s_all-elements_GSEA-like.txt" % comparison]
            tables_list.append(data)
    result = pd.concat(tables_list)
    return result


def GSEA_like_results(results_directory, tools_available):
    tool = "GSEA-like"
    if not tool in tools_available.keys():
        return "GSEA-like method wasn't activated."
    comparisons_list = os.listdir(os.path.join(results_directory, tool))
    ctrs = list(map(lambda x: x.split("_vs_")[1], list(tools_available[tool].keys())))
    @interact(control=widgets.Dropdown(options=set(ctrs), description='Control:', disabled=False),
              fdr_cutoff=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05),
              gene=widgets.Text(value='', placeholder='Feature to show...', description='Feature:'))
    def GSEA_like_results(control, fdr_cutoff, gene):
        tool = "GSEA-like"
        result = GSEA_like_data(comparison = "", control = control, tool = tool, results_directory=results_directory, tools_available=tools_available)
        result.loc[result['padj'] < fdr_cutoff, 'significant'] = 'Yes'
        result.loc[result['padj'] >= fdr_cutoff, 'significant'] = 'No'
        new_row = {'pathway':gene, 'condition':control, 'significant': 'Baseline', 'pval':1, 'padj':1 ,'ES': 0, 'NES': 0, 'size':None}
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.pathway == gene]
        domain = ['Yes', 'No', 'Baseline']
        range_ = ['red', 'grey', 'black']


        def on_button_clicked(b):
            sort_cols = natural_sort(list(res.condition.values))
            plot = alt.Chart(res).mark_circle(size=60).mark_point(
                filled=True,
                size=100,
                ).encode(
                        y='NES',
                        x=alt.X('condition:N', sort=sort_cols),
                        color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
                        tooltip=['pathway', 'condition', 'significant', 'padj', 'NES'],
                ).properties(
                        title=gene + " locfdr versus baseline (GSEA-like method)",
                        width=100
                )
            with output:
                display(plot)

        button = widgets.Button(description="Show plot")
        output = widgets.Output()
        display(button, output)
        button.on_click(on_button_clicked)



def in_house_method_data(comparison = "", control = "", tool = "", results_directory = "", tools_available = ""):
    tables_list = []
    for _comparison in os.listdir(os.path.join(results_directory, tool)):
        if _comparison.split("_vs_")[-1] == control:
            keys_list = list(tools_available[tool][_comparison].keys())
            for key in keys_list:
                if key.endswith("_all-elements_in-house.txt"):
                    break
            data = tools_available[tool][_comparison][key]
            trt = _comparison.split("_vs_")[0]
            data['condition'] = trt
            tables_list.append(data)
        if _comparison == comparison:
            data = tools_available[tool][_comparison]["%s_all-elements_in-house.txt" % comparison]
            tables_list.append(data)
    result = pd.concat(tables_list)
    return result


def in_house_method_results(results_directory, tools_available):
    tool = "in_house_method"
    if not tool in tools_available.keys():
        return "In-house method wasn't activated."
    comparisons_list = os.listdir(os.path.join(results_directory, tool))
    ctrs = list(map(lambda x: x.split("_vs_")[1], list(tools_available[tool].keys())))
    @interact(control=widgets.Dropdown(options=set(ctrs), description='Control:', disabled=False),
              fdr_cutoff=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05),
              gene=widgets.Text(value='', placeholder='Feature to show...', description='Feature:'))
    def in_house_method_results(control, fdr_cutoff, gene):
        tool = "in_house_method"
        result = in_house_method_data(comparison = "", control = control, tool = tool, results_directory=results_directory, tools_available=tools_available)
        result.loc[result['category'] == "down", 'significant'] = 'Yes'
        result.loc[result['category'] != "down", 'significant'] = 'No'
        new_row = {'Gene':gene, 'condition':control, 'significant': 'Baseline', 'down':0 ,'score': 0, 'prop': 0, 'up': 0, 'n': 0, 'category': 'Baseline'}
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.Gene == gene]
        domain = ['Yes', 'No', 'Baseline']
        range_ = ['red', 'grey', 'black']


        def on_button_clicked(b):
            sort_cols = natural_sort(list(res.condition.values))
            plot = alt.Chart(res).mark_circle(size=60).mark_point(
                filled=True,
                size=100,
                ).encode(
                        y='score',
                        x=alt.X('condition:N', sort=sort_cols),
                        color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
                        tooltip=['Gene', 'condition', 'down', 'up', 'n', 'score', 'prop'],
                ).properties(
                        title=gene + " locfdr versus baseline (Filtering method)",
                        width=100
                )
            with output:
                display(plot)

        button = widgets.Button(description="Show plot")
        output = widgets.Output()
        display(button, output)
        button.on_click(on_button_clicked)



def MAGeCK_RRA_data(comparison = "", control = "", tool = "MAGeCK_RRA", results_directory = "", tools_available = ""):

    def check(comparison, _comparison, control, mode):
        if mode:
            if _comparison.split("_vs_")[-1] == control:
                return True
            else: return False
        else:
            if _comparison == comparison:
                return True
            else: return False


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
                data['condition'] = trt
                tables_list.append(data)
    result = pd.concat(tables_list)
    return result



def MAGeCK_RRA_results(results_directory, tools_available):
    tool = "MAGeCK_RRA"
    if not tool in tools_available.keys():
        return "MAGeCK method wasn't activated."
    comparisons_list = os.listdir(os.path.join(results_directory, tool))
    ctrs = list(map(lambda x: x.split("_vs_")[1], list(tools_available[tool].keys())))
    @interact(control=widgets.Dropdown(options=set(ctrs), description='Control:', disabled=False),
              fdr_cutoff=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05),
              gene=widgets.Text(value='', placeholder='Feature to show...', description='Feature:'))
    def MAGeCK_RRA(control, fdr_cutoff, gene):
        tool = "MAGeCK_RRA"
        result = MAGeCK_RRA_data(comparison = "", control = control, tool = tool, results_directory=results_directory, tools_available=tools_available)
        result.loc[result['neg|fdr'] < fdr_cutoff, 'significant'] = 'Yes'
        result.loc[result['neg|fdr'] >= fdr_cutoff, 'significant'] = 'No'
        new_row = {'id':gene, 'condition':control, 'significant': 'Baseline', 'neg|fdr':1, 'neg|lfc': 0}
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.id == gene]
        domain = ['Yes', 'No', 'Baseline']
        range_ = ['red', 'grey', 'black']


        def on_button_clicked(b):
            sort_cols = natural_sort(list(res.condition.values))
            plot = alt.Chart(res).mark_circle(size=60).mark_point(
                filled=True,
                size=100,
                ).encode(
                        y='neg|lfc',
                        x=alt.X('condition:N', sort=sort_cols),
                        color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
                        tooltip=['id', 'neg|lfc', 'significant', 'neg|fdr', 'condition'],
                ).properties(
                        title=gene + " LFC versus baseline (MAGeCK RRA)",
                        width=100
                )
            with output:
                display(plot)

        button = widgets.Button(description="Show plot")
        output = widgets.Output()

        display(button, output)

        button.on_click(on_button_clicked)


def MAGeCK_MLE_data(comparison = "", control = "", tool = "", results_directory = "", tools_available = ""):

    def check(comparison, _comparison, control, mode):
        if mode:
            if _comparison.split("_vs_")[-1] == control:
                return True
            else: return False
        else:
            if _comparison == comparison:
                return True
            else: return False


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
                data['condition'] = trt
                if mode:
                    data = data[data.columns.drop(list(data.filter(regex= '%s' % control)))]
                    data = data.rename(columns=lambda x: re.sub('.+\|','',x))
                    tables_list.append(data)
                else:
                    tables_list.append(data)
    result = pd.concat(tables_list)
    return result


def MAGeCK_MLE_results(results_directory, tools_available):
    tool = "MAGeCK_MLE"
    if not tool in tools_available.keys():
        return "MAGeCK MLE method wasn't activated."
    comparisons_list = os.listdir(os.path.join(results_directory, tool))
    ctrs = list(map(lambda x: x.split("_vs_")[1], list(tools_available[tool].keys())))
    @interact(control=widgets.Dropdown(options=set(ctrs), description='Control:', disabled=False),
              fdr_cutoff=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05),
              gene=widgets.Text(value='', placeholder='Feature to show...', description='Feature:'))
    def MAGeCK_MLE(control, fdr_cutoff, gene):
        tool = "MAGeCK_MLE"
        result = MAGeCK_MLE_data(comparison = "", control = control, tool = tool, results_directory=results_directory, tools_available=tools_available)
        result.loc[result['fdr'] < fdr_cutoff, 'significant'] = 'Yes'
        result.loc[result['fdr'] >= fdr_cutoff, 'significant'] = 'No'
        new_row = {'Gene':gene, 'condition':control, 'significant': 'Baseline', 'fdr':1, 'beta': 0}
        result = result.append(new_row, ignore_index=True)
        res = result.loc[result.Gene == gene]

        domain = ['Yes', 'No', 'Baseline']
        range_ = ['red', 'grey', 'black']


        def on_button_clicked(b):
            sort_cols = natural_sort(list(res.condition.values))
            plot = alt.Chart(res).mark_circle(size=60).mark_point(
                filled=True,
                size=100,
                ).encode(
                        y='beta',
                        x=alt.X('condition:N', sort=sort_cols),
                        color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
                        tooltip=['Gene', 'beta', 'significant', 'fdr', 'condition'],
                ).properties(
                        title=gene + " locfdr versus baseline (MAGeCK MLE)",
                        width=100
                )
            with output:
                display(plot)

        button = widgets.Button(description="Show plot")
        output = widgets.Output()

        display(button, output)

        button.on_click(on_button_clicked)



def BAGEL_data(comparison = "", control = "", tool = "BAGEL", results_directory = "", tools_available = ""):

    def check(comparison, _comparison, control, mode):
        if mode:
            if _comparison.split("_vs_")[-1] == control:
                return True
            else: return False
        else:
            if _comparison == comparison:
                return True
            else: return False


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
                data['condition'] = trt
                tables_list.append(data)
            if _comparison == comparison:
                data = tools_available[tool][_comparison]["%s_BAGEL_output.bf" % _comparison]
                tables_list.append(data)
    result = pd.concat(tables_list)
    return result



def BAGEL_results(results_directory, tools_available):
    tool = "BAGEL"
    if not tool in tools_available.keys():
        print("BAGEL is not activated.")
    else:
        comparisons_list = os.listdir(os.path.join(results_directory, tool))
        ctrs = list(map(lambda x: x.split("_vs_")[1], list(tools_available[tool].keys())))
        @interact(control=widgets.Dropdown(options=set(ctrs), description='Control:', disabled=False),
                  fdr_cutoff=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05),
                  gene=widgets.Text(value='', placeholder='Feature to show...', description='Feature:'))
        def BAGEL(control, fdr_cutoff, gene):
            tool = "BAGEL"
            result = BAGEL_data(comparison = "", control = control, tool = tool, results_directory=results_directory, tools_available=tools_available)
            result.loc[result['BF'] > 0, 'significant'] = 'Yes'
            result.loc[result['BF'] <= 0, 'significant'] = 'No'
            new_row = {'GENE':gene, 'condition':control, 'significant': 'Baseline', 'BF':0,}
            result = result.append(new_row, ignore_index=True)
            res = result.loc[result.GENE == gene]
            domain = ['Yes', 'No', 'Baseline']
            range_ = ['red', 'grey', 'black']


            def on_button_clicked(b):
                sort_cols = natural_sort(list(res.condition.values))
                plot = alt.Chart(res).mark_circle(size=60).mark_point(
                    filled=True,
                    size=100,
                    ).encode(
                            y='BF',
                            x=alt.X('condition:N', sort=sort_cols),
                            color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
                            tooltip=["GENE", "BF", "STD", "NumObs", "condition"],
                    ).properties(
                            title=gene + " Bayes Factor versus baseline (BAGEL)",
                            width=100
                    )
                with output:
                    display(plot)

            button = widgets.Button(description="Show plot")
            output = widgets.Output()

            display(button, output)

            button.on_click(on_button_clicked)



def MAGeCK_MLE_snake_plot(comparison, fdr_cutoff, non_sig, sig, results_directory, tools_available):
    tool = "MAGeCK_MLE"
    print(tool)
    print(comparison)
    print(fdr_cutoff)
    treatment, control = comparison.split("_vs_")
    source = MAGeCK_MLE_data(comparison = comparison, control = "", tool = tool, results_directory=results_directory, tools_available=tools_available)

    source['default_rank'] = source[treatment + '|beta'].rank()
    source.loc[source[treatment + '|fdr'] < 0.05, 'significant'] = 'Yes'
    source.loc[source[treatment + '|fdr'] >= 0.05, 'significant'] = 'No'
    domain = ['Yes', 'No']
    range_ = [sig, non_sig]

    def show_chart(chart):
        try:
            print('Displaying chart...')
            display(chart)
        except:
            print("Chart can't be displayed.")

    
    def on_button_clicked(b):
        
        global chart
        
        chart = alt.Chart(source).mark_circle(size=60).encode(
            x=alt.X('default_rank:Q', axis=alt.Axis(title='Rank')),
            y=alt.Y(treatment + '|beta:Q'),
            tooltip=['Gene', 'sgRNA', treatment + '|beta', treatment + '|fdr', 'significant', 'default_rank'],
            color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
            order=alt.Order('significant:N')
        ).properties(width=800, height=400).interactive()

        line = alt.Chart(pd.DataFrame({'y': [0]})).mark_rule().encode(y='y')

        chart = (chart + line)
        show_chart(chart)
#         with output:
#             display(chart)

    button = widgets.Button(description="Show plot")
#     output = widgets.Output()
    display(button)
#     display(button, output)

    button.on_click(on_button_clicked)




def MAGeCK_RRA_snake_plot(comparison, fdr_cutoff, non_sig, sig, results_directory, tools_available):
    tool = "MAGeCK_RRA"
    print(tool)
    print(comparison)
    print(fdr_cutoff)
    treatment, control = comparison.split("_vs_")
    source = MAGeCK_RRA_data(comparison = comparison, control = "", tool = tool, results_directory=results_directory, tools_available=tools_available)

    # source['logarithm_base2'] = np.log2(source['neg|lfc']) * np.sign(source['neg|lfc'])
    source['default_rank'] = source['neg|lfc'].rank(ascending=False)
    source.loc[source['neg|fdr'] < 0.05, 'significant'] = 'Yes'
    source.loc[source['neg|fdr'] >= 0.05, 'significant'] = 'No'
    domain = ['Yes', 'No']
    range_ = [sig, non_sig]

    def on_button_clicked(b):
        chart = alt.Chart(source).mark_circle(size=60).encode(
            x=alt.X('default_rank:Q', axis=alt.Axis(title='Rank')),
            y=alt.Y('neg|lfc:Q'),
            tooltip=['id', 'num', 'neg|lfc', 'neg|fdr', 'significant', 'neg|rank', 'default_rank'],
            color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
            order=alt.Order('significant:N')
        ).properties(width=800, height=400).interactive()

        line = alt.Chart(pd.DataFrame({'y': [0]})).mark_rule().encode(y='y')

        chart = (chart + line)
        with output:
            display(chart)

    button = widgets.Button(description="Show plot")
    output = widgets.Output()

    display(button, output)

    button.on_click(on_button_clicked)



def CRISPhieRmix_snake_plot(comparison, fdr_cutoff, non_sig, sig,results_directory, tools_available):
    tool = "CRISPhieRmix"
    print(tool)
    print(comparison)
    print(fdr_cutoff)
    treatment, control = comparison.split("_vs_")
    source = CRISPhieRmix_data(comparison = comparison, control = "", tool = tool, results_directory=results_directory, tools_available=tools_available)
    source['default_rank'] = source['locfdr'].rank(method='dense')
    source.loc[source['locfdr'] < 0.05, 'significant'] = 'Yes'
    source.loc[source['locfdr'] >= 0.05, 'significant'] = 'No'
    domain = ['Yes', 'No']
    range_ = [sig, non_sig]

    def on_button_clicked(b):
        chart = alt.Chart(source).mark_circle(size=60).encode(
            x=alt.X('default_rank:Q', axis=alt.Axis(title='Rank')),
            y=alt.Y('locfdr:Q'),
            tooltip=['gene', 'locfdr', 'score', 'FDR', 'significant', 'default_rank'],
            color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
            order=alt.Order('significant:N')
        ).properties(width=800, height=400).interactive()

        line = alt.Chart(pd.DataFrame({'y': [0]})).mark_rule().encode(y='y')

        chart = (chart + line)
        with output:
            display(chart)

    button = widgets.Button(description="Show plot")
    output = widgets.Output()

    display(button, output)

    button.on_click(on_button_clicked)




def in_house_snake_plot(comparison, fdr_cutoff, non_sig, sig, results_directory, tools_available):
    tool = "in_house_method"
    print(tool)
    print(comparison)
    print(fdr_cutoff)
    treatment, control = comparison.split("_vs_")
    source = in_house_method_data(comparison = comparison, control = "", tool = tool, results_directory=results_directory, tools_available=tools_available)
    source['default_rank'] = source['score'].rank(method='dense')
    source.loc[abs(source['score']) > 1, 'significant'] = 'Yes'
    source.loc[abs(source['score']) <= 1, 'significant'] = 'No'
    domain = ['Yes', 'No']
    range_ = [sig, non_sig]

    def on_button_clicked(b):
        chart = alt.Chart(source).mark_circle(size=60).encode(
            x=alt.X('default_rank:Q', axis=alt.Axis(title='Rank')),
            y=alt.Y('score:Q'),
            tooltip=['Gene', 'up', 'down', 'n', 'significant', 'score', 'default_rank'],
            color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
            order=alt.Order('significant:N')
        ).properties(width=800, height=400).interactive()

        line = alt.Chart(pd.DataFrame({'y': [0]})).mark_rule().encode(y='y')

        chart = (chart + line)
        with output:
            display(chart)

    button = widgets.Button(description="Show plot")
    output = widgets.Output()

    display(button, output)

    button.on_click(on_button_clicked)




def GSEA_like_snake_plot(comparison, fdr_cutoff, non_sig, sig, results_directory, tools_available):
    tool = "GSEA-like"
    print(tool)
    print(comparison)
    print(fdr_cutoff)
    treatment, control = comparison.split("_vs_")
    source = GSEA_like_data(comparison = comparison, control = "", tool = tool, results_directory=results_directory, tools_available=tools_available)
    source['default_rank'] = source[['NES']].rank(method='dense')
    source.loc[abs(source['padj']) < fdr_cutoff, 'significant'] = 'Yes'
    source.loc[abs(source['padj']) >= fdr_cutoff, 'significant'] = 'No'
    domain = ['Yes', 'No']
    range_ = [sig, non_sig]

    def on_button_clicked(b):
        chart = alt.Chart(source).mark_circle(size=60).encode(
            x=alt.X('default_rank:Q', axis=alt.Axis(title='Rank')),
            y=alt.Y('NES:Q'),
            tooltip=['pathway', 'pval', 'padj', 'ES', 'NES','significant', 'size', 'default_rank'],
            color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
            order=alt.Order('significant:N')
        ).properties(width=800, height=400).interactive()

        line = alt.Chart(pd.DataFrame({'y': [0]})).mark_rule().encode(y='y')

        chart = (chart + line)
        with output:
            display(chart)

    button = widgets.Button(description="Show plot")
    output = widgets.Output()

    display(button, output)

    button.on_click(on_button_clicked)


def BAGEL_snake_plot(comparison, fdr_cutoff, non_sig, sig, results_directory, tools_available):
    tool = "BAGEL"
    print(tool)
    print(comparison)
    print(fdr_cutoff)
    treatment, control = comparison.split("_vs_")
    source = BAGEL_data(comparison = comparison, control = "", tool = tool, results_directory=results_directory, tools_available=tools_available)
    source['default_rank'] = source['BF'].rank(method='dense')
    source.loc[source['BF'] > 0, 'significant'] = 'Yes'
    source.loc[source['BF'] <= 0, 'significant'] = 'No'
    domain = ['Yes', 'No']
    range_ = [sig, non_sig]

    def on_button_clicked(b):
        chart = alt.Chart(source).mark_circle(size=60).encode(
            x=alt.X('default_rank:Q', axis=alt.Axis(title='Rank')),
            y=alt.Y('BF:Q'),
            tooltip=['GENE', 'BF', 'STD', 'significant', 'default_rank'],
            color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
            order=alt.Order('significant:N')
        ).properties(width=800, height=400).interactive()

        line = alt.Chart(pd.DataFrame({'y': [0]})).mark_rule().encode(y='y')

        chart = (chart + line)
        with output:
            display(chart)

    button = widgets.Button(description="Show plot")
    output = widgets.Output()

    display(button, output)

    button.on_click(on_button_clicked)




def snake_plot(results_directory, tools_available):
    tools = [tool for tool in tools_available.keys() if tool != "DESeq2"]
    @interact(tools=widgets.Dropdown(options=set(tools), description='Tool:', disabled=False))
    def snake_plot(tools):
        comparisons_list = os.listdir(os.path.join(results_directory, tools))
        ctrs = list(map(lambda x: x.split("_vs_")[1], list(tools_available[tools].keys())))
        @interact(comparison=widgets.Dropdown(options=comparisons_list, description='Treatment:', disabled=False),
            fdr_cutoff=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05),
            non_sig=widgets.ColorPicker(concise=False, description='Non-significant', value='#949494', disabled=False),
            sig=widgets.ColorPicker(concise=False, description='Significant', value='red', disabled=False))
        def snake_plot(comparison, fdr_cutoff, non_sig, sig):
            if tools == "MAGeCK_RRA":
                MAGeCK_RRA_snake_plot(comparison, fdr_cutoff, non_sig, sig, results_directory, tools_available)
            if tools == "MAGeCK_MLE":
                MAGeCK_MLE_snake_plot(comparison, fdr_cutoff, non_sig, sig, results_directory, tools_available)
            if tools == "CRISPhieRmix":
                CRISPhieRmix_snake_plot(comparison, fdr_cutoff, non_sig, sig, results_directory, tools_available)
            if tools == "in_house_method":
                in_house_snake_plot(comparison, fdr_cutoff, non_sig, sig, results_directory, tools_available)
            if tools == "GSEA-like":
                GSEA_like_snake_plot(comparison, fdr_cutoff, non_sig, sig, results_directory, tools_available)
            if tools == "BAGEL":
                BAGEL_snake_plot(comparison, fdr_cutoff, non_sig, sig, results_directory, tools_available)
            else:
                print("Choose a tool.")


def show_sgRNA_counts(token):
    config = "./config/%s.yaml" % token
    content = open_yaml(config)
    cts_file = content['normalized_count_table']
    cts = pd.read_csv(cts_file, sep="\t")
    cts_columns = [col for col in cts.columns.tolist() if not col in ["sgRNA", "Gene"]]
    @interact(element=widgets.Text(value='', placeholder='Element:', description='Element:'),
              conditions=widgets.TagsInput(value=cts_columns, allowed_tags=cts_columns, allow_duplicates=False))
    def show_sgRNA_counts(element, conditions):
        button = widgets.Button(description="Show sgRNA counts...")
        output = widgets.Output()
        display(button, output)
        def on_button_clicked(b):
            cts = pd.read_csv(cts_file, sep="\t")
            cts = pd.melt(cts, id_vars=['sgRNA', 'Gene'])
            if not element in list(cts.Gene):
                gene_cts = cts
            else:
                gene_cts = cts.loc[cts.Gene == element]
            gene_cts = gene_cts.loc[gene_cts.variable.isin(conditions)]
            gene_cts = gene_cts.pivot_table(index="sgRNA", columns="variable", values="value")
            gene_cts = gene_cts[conditions]
            net.load_df(gene_cts)
            net.normalize(axis='row', norm_type='zscore')
            net.cluster()
            with output:
                display(net.widget())

        button.on_click(on_button_clicked)


def show_sgRNA_counts_lines(token):
    config = "./config/%s.yaml" % token
    content = open_yaml(config)
    cts_file = content['normalized_count_table']
    cts = pd.read_csv(cts_file, sep="\t")
    design_file = content['tsv_file']
    design = pd.read_csv(design_file, sep="\t")
    @interact(element=widgets.Text(value='', placeholder='Element:', description='Element:'),
              conditions=widgets.Text(value=",".join(list(set(design.condition))), placeholder='Conditions to show, in order and comma-separated:', description='Conditions:'))
    def show_sgRNA_counts_lines(element, conditions):
        sort_cols = conditions.split(",")
        button = widgets.Button(description="Show sgRNA counts...")
        output = widgets.Output()
        display(button, output)
        def on_button_clicked(b):
            cts = pd.read_csv(cts_file, sep="\t")
            cts = pd.melt(cts, id_vars=['sgRNA', 'Gene'])
            if not element in list(cts.Gene):
                print("Element '%s' not in counts matrix." % element)
            if not len(conditions.split(",")) > 1:
                print("Choose more conditions to show.")
            gene_cts = cts.loc[cts.Gene == element]
            source = gene_cts
            source = pd.merge(source, design, left_on='variable', right_on='replicate')
            source = source.groupby(['sgRNA', 'condition']).mean()
            source = source.reset_index()
            boolean_series = source.condition.isin(sort_cols)
            source = source[boolean_series]
            out = alt.Chart(source).mark_line().encode(
                x=alt.X('condition', axis=alt.Axis(title='Condition'), sort=sort_cols),
                y='value',
                color='sgRNA'
            ).interactive().properties(width=300)
            with output:
                display(out)
        button.on_click(on_button_clicked)


### Enrichr
def getEnrichrResults(genes, description, gene_set_library):
    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(genes)
    description = description
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)

    userListId = data['userListId']

    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = userListId
    gene_set_library = gene_set_library[:-1]
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
     )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text)
    return data


def createEnrichrTable(enrichrResults):
    cols = ['Rank', 'Term name', 'P-value', 'Z-score', 'Combined score', 'Overlapping genes', 'Adjusted p-value', 'Old p-value', 'Old adjusted p-value']
    for base in enrichrResults:
        rows = []
        for pathway in enrichrResults[base]:
            row = pd.DataFrame([pathway])
            rows.append(row)

    table = pd.concat(rows)
    table.columns = cols
    return table


def enrichmentBarPlot(source, n, description, col_1, col_2, base):
    if n == 'max':
        n = len(source.index)
    source = source.sort_values(by=['Combined score'], ascending=False).head(n)

    domain = [source['Adjusted p-value'].min(), source['Adjusted p-value'].max()]
    range_ = [col_1, col_2]


    bars = alt.Chart(source).mark_bar().encode(
        x='Combined score',
        y=alt.Y('Term name', sort='-x'),
        tooltip=['Rank',
                 'Term name',
                 'P-value',
                 'Z-score',
                 'Combined score',
                 'Overlapping genes',
                 'Adjusted p-value',
                 'Old p-value',
                 'Old adjusted p-value'],
        color=alt.Color('Adjusted p-value', scale=alt.Scale(domain=domain, range=range_),
                        legend=alt.Legend(title="Adjusted p-value:")),
    ).properties(
        title=description + " (%s)" % base[:-1],
    )


    chart = (bars).properties(height=15*n, width=500)
    return chart


def enrichmentCirclePlot(source, n, description, col_1, col_2, base):
    if n == 'max':
        n = int(len(source.index))

    source = source.sort_values(by=['Combined score'], ascending=False).head(n)

    source['Overlap size'] = source['Overlapping genes'].str.len()

    domain = [source['Adjusted p-value'].min(), source['Adjusted p-value'].max()]
    range_ = [col_1, col_2]

    chart = alt.Chart(source).mark_circle(size=60).encode(
        alt.X('Combined score',
            scale=alt.Scale(domain=(source.min()['Combined score'], source.max()['Combined score']*1.05))),
        y=alt.Y('Term name', sort='-x'),
        color=alt.Color('Adjusted p-value', scale=alt.Scale(domain=domain, range=range_),
                        legend=alt.Legend(title="Adjusted p-value:")),
        tooltip=['Rank',
                 'Term name',
                 'P-value',
                 'Z-score',
                 'Combined score',
                 'Overlapping genes',
                 'Adjusted p-value',
                 'Old p-value',
                 'Old adjusted p-value',
                 'Overlap size'],
        size='Overlap size',
    ).properties(
        title=description + " (%s)" % base[:-1]
    ).configure_axisY(
        titleAngle=0,
        titleY=-10,
        titleX=-60,
        labelLimit=250
    ).interactive()

    chart = (chart).properties(height=20*n, width=500)

    return chart







def enrichr_plots(pitviper_res):
    BASES = open("workflow/notebooks/enrichr_list.txt", "r").readlines()
    TOOLS = [tool for tool in pitviper_res.keys() if tool != "DESeq2"]
    @interact(tool=TOOLS)
    def enrichr_plots(tool):
        tool_res = pitviper_res[tool]
        conditions = tool_res.keys()
        @interact(description=widgets.Text(value='My gene list', placeholder='Description', description='Description:'), bases=widgets.SelectMultiple(options=BASES), conditions=conditions)
        def enrichr_plots(description, bases, conditions):
            treatment, baseline = conditions.split("_vs_")
            @interact(col_2=widgets.ColorPicker(concise=False, description='Top color', value='blue', disabled=False), col_1=widgets.ColorPicker(concise=False, description='Bottom color', value='red', disabled=False), plot_type=['Circle', 'Bar'], size=[5, 10, 20, 50, 100, 200, 'max'], fdr_cutoff=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05), score_cutoff=widgets.Text(value='0', placeholder='0', description='Score cut-off:'))
            def enrichr_plots(score_cutoff, fdr_cutoff, size, plot_type, col_2, col_1):
                print('Description:', description)
                print('Baseline:', baseline)
                print('Treatment:', treatment)
                print('FDR cut-off:', fdr_cutoff)
                print('Score cut-off', score_cutoff)
                print('Tool:', tool)
                print('Gene set library:', bases)

                score_cutoff = float(score_cutoff)

                def on_button_clicked(b):
                    charts = []
                    for base in bases:
                        if tool == 'MAGeCK_MLE':
                            info = tool_res[conditions][conditions + ".gene_summary.txt"]
                            info = info.loc[info[treatment+'|fdr'] < fdr_cutoff]
                            if score_cutoff < 0:
                                info = info.loc[info[treatment+'|beta'] < score_cutoff]
                            elif score_cutoff > 0:
                                info = info.loc[info[treatment+'|beta'] > score_cutoff]
                            genes = info['Gene']

                        if tool == 'MAGeCK_RRA':
                            info = tool_res[conditions][conditions + ".gene_summary.txt"]
                            info = info.loc[info['neg|fdr'] < fdr_cutoff]
                            genes = info['id']

                        if tool == 'BAGEL':
                            info = tool_res[conditions][conditions + "_BAGEL_output.bf"]
                            info = info.loc[info['BF'] > score_cutoff]
                            genes = info['GENE']

                        if tool == 'in_house_method':
                            info = tool_res[conditions][conditions + "_all-elements_in-house.txt"]
                            info = info.loc[info['score'] < score_cutoff]
                            genes = info['Gene']

                        if tool == "GSEA-like":
                            info = tool_res[conditions][conditions + "_all-elements_GSEA-like.txt"]
                            if score_cutoff > 0:
                                info = info.loc[info['NES'] > score_cutoff]
                            elif score_cutoff < 0:
                                info = info.loc[info['NES'] < score_cutoff]
                            info = info.loc[info['padj'] < fdr_cutoff]
                            genes = info['pathway']

                        # print("Size (gene set):", len(genes))

                        enrichr_res = getEnrichrResults(genes, description, base)
                        table = createEnrichrTable(enrichr_res)
                        if plot_type == 'Bar':
                            chart = enrichmentBarPlot(table, size, description, col_1, col_2, base)
                        else:
                            chart = enrichmentCirclePlot(table, size, description, col_1, col_2, base)
                        charts.append(chart)
                    with output:
                        for chart in charts:
                            display(chart)


                button = widgets.Button(description="Show EnrichR results")
                output = widgets.Output()

                display(button, output)

                button.on_click(on_button_clicked)


def pca_counts(token):
    config = "./config/%s.yaml" % token
    content = open_yaml(config)

    TSV = pd.read_csv(content['tsv_file'], sep="\t")
    cts_file = content['normalized_count_table']
    cts = pd.read_csv(cts_file, sep="\t")
    X = cts[cts.columns[2:]].to_numpy().T
    d = dict(zip(TSV.replicate, TSV.condition))
    y = [d[k] for k in cts.columns[2:]]
    y = np.array(y)
    y_bis = np.array(cts.columns[2:])

    pca = decomposition.PCA(n_components=2)
    pca.fit(X)
    X = pca.transform(X)

    a = pd.DataFrame(X, columns=['dim1', 'dim2'])
    b = pd.DataFrame(y, columns=['condition'])
    c = pd.DataFrame(y_bis, columns=['replicate'])

    df_c = pd.concat([a, b, c], axis=1)

    source = df_c

    pca_2d = alt.Chart(source).mark_circle(size=60).encode(
        x='dim1',
        y='dim2',
        color='condition:N',
        tooltip=['dim1', 'dim2', 'condition', 'replicate']
    ).interactive()

    return pca_2d


def ranking(treatment, control, token, tools_available, params):

    def get_occurence_df(data):
        essential_genes = []
        for key in list(data.keys()):
            essential_genes.extend(data[key])
        df = pd.DataFrame(np.zeros((len(set(essential_genes)), len(data.keys()))),
                          index = set(essential_genes),
                          columns = data.keys(),
                          dtype = int)
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

    if params['MAGeCK_MLE']['on']:
        score = params['MAGeCK_MLE']['score']
        fdr = params['MAGeCK_MLE']['fdr']
        greater = params['MAGeCK_MLE']['greater']
        mle = tools_available["MAGeCK_MLE"][comparison][comparison + ".gene_summary.txt"]
        if not greater:
            mle = mle[(mle["%s|beta" % treatment] < score) & (mle["%s|fdr" % treatment] < fdr)]
        else:
            mle = mle[(mle["%s|beta" % treatment] > score) & (mle["%s|fdr" % treatment] < fdr)]
        mle['default_rank'] = mle[treatment + '|beta'].rank(method="dense").copy()
        mle = mle[["Gene", "default_rank"]].rename(columns={"Gene": "id", "default_rank": "mle_rank"})
        mle_genes = list(mle.id)
        tool_results["MAGeCK MLE"] = mle_genes
        tool_genes.append(mle_genes)

    if params['MAGeCK_RRA']['on']:
        score = params['MAGeCK_RRA']['score']
        fdr = params['MAGeCK_RRA']['fdr']
        greater = params['MAGeCK_RRA']['greater']
        direction = params['MAGeCK_RRA']['direction']
        rra = tools_available["MAGeCK_RRA"][comparison][comparison + ".gene_summary.txt"]
        if direction == 'Negative':
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


    if params['BAGEL']['on']:
        score = params['BAGEL']['score']
        greater = params['BAGEL']['greater']
        bagel = tools_available["BAGEL"][comparison][comparison + "_BAGEL_output.bf"]
        if greater:
            bagel = bagel[(bagel["BF"] > score)]
        else:
            bagel = bagel[(bagel["BF"] < score)]
        bagel['default_rank'] = bagel['BF'].rank(method="dense", ascending=False).copy()
        bagel = bagel[["GENE", "default_rank"]].rename(columns={"GENE": "id", "default_rank": "bagel_rank"})
        bagel_genes = list(bagel.id)
        tool_results["BAGEL"] = bagel_genes
        tool_genes.append(bagel_genes)

    if params['in_house_method']['on']:
        score = params['in_house_method']['score']
        greater = params['in_house_method']['greater']
        in_house = tools_available["in_house_method"][comparison][comparison + "_all-elements_in-house.txt"]
        if params['in_house_method']['direction'] == "Negative":
            if not greater:
                in_house = in_house[(in_house["down"] > 1) & (in_house["up"] < 2) & (in_house["score"] < score)]
            else:
                in_house = in_house[(in_house["down"] > 1) & (in_house["up"] < 2) & (in_house["score"] > score)]
        elif params['in_house_method']['direction'] == "Positive":
            if not greater:
                in_house = in_house[(in_house["down"] < 2) & (in_house["up"] > 1) & (in_house["score"] < score)]
            else:
                in_house = in_house[(in_house["down"] < 2) & (in_house["up"] > 1) & (in_house["score"] > score)]
        in_house['default_rank'] = in_house['score'].rank(method="dense").copy()
        in_house = in_house[["Gene", "default_rank"]].rename(columns={"Gene": "id", "default_rank": "in_house_rank"})
        in_house_genes = list(in_house.id)
        tool_results["In-house"] = in_house_genes
        tool_genes.append(in_house_genes)


    if params['GSEA_like']['on']:
        score = params['GSEA_like']['score']
        fdr = params['GSEA_like']['fdr']
        greater = params['GSEA_like']['greater']
        gsea = tools_available["GSEA-like"][comparison][comparison + "_all-elements_GSEA-like.txt"]
        if not greater:
            gsea = gsea[(gsea["NES"] < score) & (gsea["pval"] < fdr)]
        else:
            gsea = gsea[(gsea["NES"] > score) & (gsea["pval"] < fdr)]
        gsea['default_rank'] = gsea['NES'].rank(method="dense").copy()
        gsea = gsea[["pathway", "default_rank"]].rename(columns={"pathway": "id", "default_rank": "gsea_rank"})
        gsea_genes = list(gsea.id)
        tool_results["GSEA-like"] = gsea_genes
        tool_genes.append(gsea_genes)


    if params['CRISPhieRmix']['on']:
        score = params['CRISPhieRmix']['score']
        fdr = params['CRISPhieRmix']['fdr']
        greater = params['CRISPhieRmix']['greater']
        crisphie = tools_available["CRISPhieRmix"][comparison][comparison + ".txt"]
        if not greater:
            crisphie = crisphie[(crisphie["score"] < score) & (crisphie["FDR"] < fdr)]
        else:
            crisphie = crisphie[(crisphie["score"] < score) & (crisphie["FDR"] < fdr)]
        crisphie['default_rank'] = crisphie['FDR'].rank(method="dense").copy()
        crisphie = crisphie[["gene", "default_rank"]].rename(columns={"gene": "id", "default_rank": "crisphiermix_rank"})
        crisphie_genes = list(crisphie.id)
        tool_results["CRISPhieRmix"] = crisphie_genes
        tool_genes.append(crisphie_genes)

    l = []
    for genes in tool_genes:
        for gene in genes:
            l.append(gene)

    if params['MAGeCK_MLE']['on']:
        mle = tools_available["MAGeCK_MLE"][comparison][comparison + ".gene_summary.txt"]
        mle['default_rank'] = mle[treatment + '|beta'].rank(method="dense").copy()
        mle = mle[["Gene", "default_rank"]].rename(columns={"Gene": "id", "default_rank": "mle_rank"})
        pdList.append(mle)

    if params['MAGeCK_RRA']['on']:
        rra = tools_available["MAGeCK_RRA"][comparison][comparison + ".gene_summary.txt"]
        if params['MAGeCK_RRA']['direction'] == 'Negative':
            rra = rra[["id", "neg|rank"]].rename(columns={"neg|rank": "rra_rank"})
        elif params['MAGeCK_RRA']['direction'] == 'Positive':
            rra = rra[["id", "pos|rank"]].rename(columns={"pos|rank": "rra_rank"})
        pdList.append(rra)

    if params['BAGEL']['on']:
        bagel = tools_available["BAGEL"][comparison][comparison + "_BAGEL_output.bf"]
        bagel['default_rank'] = bagel['BF'].rank(method="dense", ascending=False).copy()
        bagel = bagel[["GENE", "default_rank"]].rename(columns={"GENE": "id", "default_rank": "bagel_rank"})
        pdList.append(bagel)

    if params['in_house_method']['on']:
        in_house = tools_available["in_house_method"][comparison][comparison + "_all-elements_in-house.txt"]
        in_house['default_rank'] = in_house['score'].rank(method="dense").copy()
        in_house = in_house[["Gene", "default_rank"]].rename(columns={"Gene": "id", "default_rank": "in_house_rank"})
        pdList.append(in_house)

    if params['GSEA_like']['on']:
        gsea = tools_available["GSEA-like"][comparison][comparison + "_all-elements_GSEA-like.txt"]
        gsea['default_rank'] = gsea['NES'].rank(method="dense").copy()
        gsea = gsea[["pathway", "default_rank"]].rename(columns={"pathway": "id", "default_rank": "gsea_rank"})
        pdList.append(gsea)

    if params['CRISPhieRmix']['on']:
        crisphie = tools_available["CRISPhieRmix"][comparison][comparison + ".txt"]
        crisphie['default_rank'] = crisphie['FDR'].rank(method="dense").copy()
        crisphie = crisphie[["gene", "default_rank"]].rename(columns={"gene": "id", "default_rank": "crisphiermix_rank"})
        print(crisphie.head())
        pdList.append(crisphie)

    df_merged_reduced = reduce(lambda  left,right: pd.merge(left,right,on=['id'],
                                                how='outer'), pdList)

    df_merged_reduced = df_merged_reduced[df_merged_reduced['id'].isin(l)]

    occurence_df = get_occurence_df(tool_results)

    return (df_merged_reduced, occurence_df)


def genemania_link_results(tools_available):
    TOOLS = [tool for tool in tools_available.keys() if tool != "DESeq2"]
    @interact(tool=TOOLS)
    def genemania_link_results(tool):
        tool_res = tools_available[tool]
        conditions = tool_res.keys()
        @interact(conditions=conditions)
        def genemania_link_results(conditions):
            treatment, baseline = conditions.split("_vs_")
            @interact(fdr_cutoff=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05),
                      score_cutoff=widgets.FloatText(value=0, description='Score cut-off:'))
            def genemania_link_results( fdr_cutoff, score_cutoff):
                print('Baseline:', baseline)
                print('Treatment:', treatment)
                print('FDR cut-off:', fdr_cutoff)
                print('Score cut-off', score_cutoff)
                print('Tool:', tool)
                html = ""

                def on_button_clicked(b):
                    if tool == 'MAGeCK_MLE':
                        info = tool_res[conditions][conditions + ".gene_summary.txt"]
                        info = info.loc[info[treatment+'|fdr'] < fdr_cutoff]
                        if float(score_cutoff) < 0:
                            info = info.loc[info[treatment+'|beta'] < score_cutoff]
                        elif float(score_cutoff) > 0:
                            info = info.loc[info[treatment+'|beta'] > score_cutoff]
                        genes = info['Gene']

                    if tool == 'MAGeCK_RRA':
                        info = tool_res[conditions][conditions + ".gene_summary.txt"]
                        info = info.loc[info['neg|fdr'] < fdr_cutoff]
                        genes = info['id']

                    if tool == 'BAGEL':
                        info = tool_res[conditions][conditions + "_BAGEL_output.bf"]
                        info = info.loc[info['BF'] > float(score_cutoff)]
                        genes = info['GENE']

                    if tool == 'in_house_method':
                        info = tool_res[conditions][conditions + "_all-elements_in-house.txt"]
                        info = info.loc[info['score'] < float(score_cutoff)]
                        genes = info['Gene']

                    if tool == "GSEA-like":
                        info = tool_res[conditions][conditions + "_all-elements_GSEA-like.txt"]
                        if float(score_cutoff) > 0:
                            info = info.loc[info['NES'] > float(score_cutoff)]
                        elif float(score_cutoff) < 0:
                            info = info.loc[info['NES'] < float(score_cutoff)]
                        info = info.loc[info['padj'] < fdr_cutoff]
                        genes = info['pathway']

                    with output:
                        print("Size (gene set):", len(genes))
                        link = "http://genemania.org/search/homo-sapiens/" + "/".join(genes)
                        print(link)

                button = widgets.Button(description="Genemania")
                output = widgets.Output()

                display(button, output)

                button.on_click(on_button_clicked)


def plot_venn(occurences):
    venn_lib = importr('venn')
    grdevices = importr('grDevices')
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_occurences = ro.conversion.py2rpy(occurences)

    with rpy2.robjects.lib.grdevices.render_to_bytesio(grdevices.png, width=1024, height=896, res=150) as img:
        venn_lib.venn(r_occurences, ilabels = False, zcolor = "style", ilcs= 1, sncs = 1, borders = False, box = False)
    IPython.display.display(IPython.display.Image(data=img.getvalue(), format='png', embed=True, height=896, width=512))


def run_rra(ranks):
    rra_lib = importr('RobustRankAggreg')
    r_source = robjects.r['source']
    r_source("workflow/notebooks/functions_R.R")
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_ranks = ro.conversion.py2rpy(ranks)
    RobustRankAggregate = robjects.r['RobustRankAggregate']
    res = RobustRankAggregate(r_ranks)
    print(res)


def reset_params():
    params = {'MAGeCK_MLE': {'on': False,
                             'fdr': 0.05,
                             'score': 0,
                             'greater': True},
              'MAGeCK_RRA': {'on': False,
                             'fdr': 0.05,
                             'score': 0,
                             'greater': True,
                             'direction': 'Positive'},
              'BAGEL': {'on': False,
                             'score': 0,
                             'greater': True},
              'CRISPhieRmix': {'on': False,
                             'fdr': 0.05,
                             'score': 0,
                             'greater': True},
              'in_house_method': {'on': False,
                             'direction': 'Positive',
                             'score': 0,
                             'greater': True},
              'GSEA_like': {'on': False,
                             'fdr': 0.05,
                             'score': 0,
                             'greater': True}}
    return params



def intersection(tools_available, token):
    TOOLS = [tool for tool in tools_available.keys() if not tool in ["DESeq2", "CRISPhieRmix"]]

    # Define widgets's options
    conditions_options = tools_available[TOOLS[0]].keys()
    tools_options = TOOLS

    # Define widgets
    conditions_widget = widgets.Dropdown(options=conditions_options)
    tools_widget = widgets.SelectMultiple(options=tools_options)

    # Show form selection
    def run(condition_selected, tools_selected):

        params = reset_params()
        global enrichr_tools
        enrichr_tools = ()

        def mle_order_update(change):
            with output:
                if change['new'] == 'Greater than score':
                    params['MAGeCK_MLE']['greater'] = True
                else:
                    params['MAGeCK_MLE']['greater'] = False
        def mle_fdr_update(change):
            with output:
                params['MAGeCK_MLE']['fdr'] = change['new']
        def mle_score_update(change):
            with output:
                params['MAGeCK_MLE']['score'] = change['new']

        def rra_order_update(change):
            with output:
                if change['new'] == 'Greater than score':
                    params['MAGeCK_RRA']['greater'] = True
                else:
                    params['MAGeCK_RRA']['greater'] = False
        def rra_fdr_update(change):
            with output:
                params['MAGeCK_RRA']['fdr'] = change['new']
        def rra_score_update(change):
            with output:
                params['MAGeCK_RRA']['score'] = change['new']
        def rra_direction_update(change):
            with output:
                params['MAGeCK_RRA']['direction'] = change['new']

        def bagel_order_update(change):
            with output:
                if change['new'] == 'Greater than score':
                    params['BAGEL']['greater'] = True
                else:
                    params['BAGEL']['greater'] = False
        def bagel_fdr_update(change):
            with output:
                params['BAGEL']['fdr'] = change['new']
        def bagel_score_update(change):
            with output:
                params['BAGEL']['score'] = change['new']

        def CRISPhieRmix_order_update(change):
            with output:
                if change['new'] == 'Greater than score':
                    params['CRISPhieRmix']['greater'] = True
                else:
                    params['CRISPhieRmix']['greater'] = False
        def CRISPhieRmix_fdr_update(change):
            with output:
                params['CRISPhieRmix']['fdr'] = change['new']
        def CRISPhieRmix_score_update(change):
            with output:
                params['CRISPhieRmix']['score'] = change['new']

        def in_house_method_order_update(change):
            with output:
                if change['new'] == 'Greater than score':
                    params['in_house_method']['greater'] = True
                else:
                    params['in_house_method']['greater'] = False
        def in_house_method_direction_update(change):
            with output:
                params['in_house_method']['direction'] = change['new']
        def in_house_method_score_update(change):
            with output:
                params['in_house_method']['score'] = change['new']


        def GSEA_like_order_update(change):
            with output:
                if change['new'] == 'Greater than score':
                    params['GSEA_like']['greater'] = True
                else:
                    params['GSEA_like']['greater'] = False
        def GSEA_like_fdr_update(change):
            with output:
                params['GSEA_like']['fdr'] = change['new']
        def GSEA_like_score_update(change):
            with output:
                params['GSEA_like']['score'] = change['new']

        def venn_button_clicked(b):
            with output:
                clear_output(wait=True)
                global treatment
                global control
                treatment, control = condition_selected.split("_vs_")
                global ranks
                global occurences
                ranks, occurences = ranking(treatment, control, token, tools_available, params)
                df = pd.DataFrame(occurences.eq(occurences.iloc[:, 0], axis=0).all(1), columns = ['intersection'])
                intersection_genes = df.loc[df.intersection == True].index
                plot_venn(occurences)
                print("Genes at intersection of all methods:")
                for gene in intersection_genes:
                    print(gene)


        def rra_button_clicked(b):
            with output:
                clear_output(wait=True)
                global treatment
                global control
                treatment, control = condition_selected.split("_vs_")
                global ranks
                global occurences
                ranks, occurences = ranking(treatment, control, token, tools_available, params)
                run_rra(ranks)


        def genemania_button_clicked(b):
            with output:
                clear_output(wait=True)
                global treatment
                global control
                treatment, control = condition_selected.split("_vs_")
                global ranks
                global occurences
                ranks, occurences = ranking(treatment, control, token, tools_available, params)
                genes_list = list(occurences[occurences.all(axis='columns')].index)
                link = "http://genemania.org/search/homo-sapiens/" + "/".join(genes_list)
                print('Link to Genemania website: (%s elements)\n' % len(genes_list))
                print(link)



        def enrichr_button_clicked(b):
            def test(b, genes, bases, size, plot_type, col_2, col_1, description):
                with output:
                    charts = []
                    for base in bases:
                        enrichr_res = getEnrichrResults(genes, description, base)
                        table = createEnrichrTable(enrichr_res)
                        if plot_type == 'Bar':
                            chart = enrichmentBarPlot(table, size, description, col_1, col_2, base)
                        else:
                            chart = enrichmentCirclePlot(table, size, description, col_1, col_2, base)
                        charts.append(chart)
                    for chart in charts:
                        display(chart)
            with output:
                clear_output(wait=True)
                treatment, control = condition_selected.split("_vs_")
                ranks, occurences = ranking(treatment, control, token, tools_available, params)
                genes_list = list(occurences[occurences.all(axis='columns')].index)
                BASES = open("workflow/notebooks/enrichr_list.txt", "r").readlines()
                @interact(bases = widgets.SelectMultiple(options=BASES, description='Bases:'),
                          size = widgets.Dropdown(options = [5, 10, 20, 50, 100, 200, 'max'], description='Size:'),
                          plot_type = widgets.Dropdown(options = ['Circle', 'Bar'], description='Type:'),
                          col_2 = widgets.ColorPicker(concise=False, description='Top color', value='blue', disabled=False),
                          col_1 = widgets.ColorPicker(concise=False, description='Bottom color', value='red', disabled=False),
                          description = widgets.Text(value='My gene list', placeholder='Description', description='Description:'))
                def enrichr_interact(bases, size, plot_type, col_2, col_1, description):
                    button_enrichr = widgets.Button(description="Run EnrichR!")
                    display(button_enrichr)
                    button_enrichr.on_click(functools.partial(test, genes=genes_list, bases=bases, size=size, plot_type=plot_type, col_2=col_2, col_1=col_1, description=description))




        def depmap_button_clicked(b):
            data_types_widget = widgets.RadioButtons(options=['crispr', 'proteomic', 'rnai', 'tpm', 'mutations'],value='crispr',description='Data type:',disabled=False)

            
            def download_depmap_file(data_type, release):
                target_file = "resources\/depmap\/%s_%s.txt" % (release, data_type)
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

                    
            def getRelease():
                depmap_release = depmap.depmap_release()
#                 utils.write_table(depmap_release, "resources/depmap/release.txt", row_names=False, quote=False, col_names=False)
#                 release = open("resources/depmap/release.txt", "r").read()
                return str(release).rstrip()
                
    
            def depmap_query_button_clicked(b):
                depmap_release = getRelease()
                save_path = "resources/depmap/%s_%s.txt" % (depmap_release, data_types_widget.value)
                Path("resources/depmap/").mkdir(parents=True, exist_ok=True)
                treatment, control = condition_selected.split("_vs_")
                ranks, occurences = ranking(treatment, control, token, tools_available, params)
                with output:
                    clear_output(wait=True)
                    if download_depmap_file(data_types_widget.value, depmap_release):
                        print("This step can take some time.")
                        print("Querying: %s..." % data_types_widget.value)
                        eh = experimentHub.ExperimentHub()
                        base_package = importr("base")
                        dplyr = importr("dplyr")
                        tidyr = importr("tidyr")
                        readr = importr("readr")
                        if data_types_widget.value == "rnai":
                            depmap_data = tidyr.drop_na(dplyr.select(depmap.depmap_rnai(), "depmap_id", "gene_name", "dependency"))
                        elif data_types_widget.value == "crispr":
                            depmap_data = tidyr.drop_na(dplyr.select(depmap.depmap_crispr(), "depmap_id", "gene_name", "dependency"))
                        elif data_types_widget.value == "proteomic":
                            depmap_data = tidyr.drop_na(dplyr.select(depmap.depmap_proteomic(), "depmap_id", "gene_name", "protein_expression"))
                        elif data_types_widget.value == "tpm":
                            depmap_data = tidyr.drop_na(dplyr.select(depmap.depmap_TPM(), "depmap_id", "gene_name", "rna_expression"))
                        elif data_types_widget.value == "mutations":
                            depmap_data = tidyr.drop_na(dplyr.select(depmap.depmap_mutationCalls(), 'depmap_id', 'gene_name', 'protein_change', 'is_deleterious'))
                        if not os.path.isfile("resources/depmap/%s_metadata.txt" % depmap_release):
                            depmap_metadata = dplyr.select(depmap.depmap_metadata(), 'depmap_id', 'sample_collection_site', 'primary_or_metastasis', 'primary_disease','subtype_disease', 'cell_line', 'cell_line_name')
                            print("Saving metadata...")
                            utils.write_table(depmap_metadata, "resources/depmap/%s_metadata.txt" % depmap_release, row_names=False, quote=False, sep="\t")
                        else:
                            print("Import metadata...")
                            depmap_metadata = readr.read_delim("resources/depmap/%s_metadata.txt" % depmap_release, delim="\t")
                        depmap_data = base_package.merge(depmap_data, depmap_metadata, by = "depmap_id")
                        print("Saving %s" % save_path)
                        utils.write_table(depmap_data, save_path, row_names=False, quote=False, sep="\t")
                    print("Opening %s" % save_path)
                    py_tissues = pd.read_table(save_path, sep="\t")
                    cell_line = list(set(py_tissues.cell_line))
                    tissues = ["_".join(str(tissu).split("_")[1:]) for tissu in cell_line]
                    tissues = list(set(tissues))
                    tissues.insert(0, 'All')
                    tissues_widget = widgets.SelectMultiple(options = tissues, description="Tissu:", value=(tissues[0], ))
                    def tissu_selection_button_clicked(b):
                        dic = {"rnai":"dependency", "crispr":"dependency", "tpm":"rna_expression", "proteomic":"protein_expression", "mutations":["protein_change","is_deleterious"]}
                        variable = dic[data_types_widget.value]
                        save_path = "resources/depmap/%s_%s.txt" % (depmap_release, data_types_widget.value)
                        py_tissues = pd.read_table(save_path, sep="\t")
                        table = py_tissues[py_tissues['primary_disease'].str.contains('|'.join(primary_tissu_widget.value), na=False)]
                        table = table[table['cell_line'].str.contains('|'.join(cell_lines_widget.value), na=False)]
                        genes_list = list(occurences[occurences.all(axis='columns')].index)
                        if data_types_widget.value == "mutations":
                            columns = ['gene_name', 'cell_line', 'cell_line_name', 'sample_collection_site', 'primary_or_metastasis', 'primary_disease','subtype_disease']
                            for value in variable:
                                columns.append(value)
                            essential_genes = table[table.gene_name.isin(genes_list)][columns]
                            essential_genes = essential_genes.loc[essential_genes.is_deleterious == True]
                            essential_genes = essential_genes.groupby(['gene_name', 'cell_line', 'cell_line_name', 'sample_collection_site', 'primary_or_metastasis', 'primary_disease',
       'subtype_disease'])['protein_change'].apply(lambda x: "%s" % '|'.join(map(str, x))).reset_index()
                            chart = alt.Chart(
                                essential_genes,
                                title="depmap mutations"
                            ).mark_rect().encode(
                                x='cell_line_name',
                                y='gene_name',
                                color=alt.Color('primary_disease', scale=alt.Scale(scheme="tableau20")),
                                tooltip=['protein_change', 'gene_name', 'cell_line', 'cell_line_name', 'sample_collection_site', 'primary_or_metastasis', 'primary_disease','subtype_disease']
                            ).interactive()
                            display(chart)
                        else:
                            essential_genes = table[table.gene_name.isin(genes_list)][["gene_name", "cell_line", variable]].pivot_table(index='gene_name', columns='cell_line', values=variable).dropna()
                            essential_genes = essential_genes.rename(columns=lambda s: s.split("_")[0])
                            net.load_df(essential_genes)
                            if data_types_widget.value == "tpm":
                                net.normalize(axis='row', norm_type='zscore')
                            net.cluster()
                            display(net.widget())
                                
                    if tissues_widget.value != "All":
                        primary_tissu_list = list(set(py_tissues[py_tissues['cell_line'].str.contains('|'.join(tissues_widget.value), na=False)].primary_disease))
                    else:
                        primary_tissu_list = list(set(py_tissues['primary_disease'].tolist()))
                    primary_tissu_list.insert(0, 'All')
                    primary_tissu_widget = widgets.SelectMultiple(options = primary_tissu_list, description="Primary tissu:", value=(primary_tissu_list[0],))
                    if primary_tissu_widget.value != "All":
                        cell_lines_list = list(set(py_tissues[py_tissues['primary_disease'].str.contains('|'.join(primary_tissu_widget.value), na=False)].cell_line_name))
                    else:
                        cell_lines_list = list(set(py_tissues['cell_line_name'].to_list()))
                    cell_lines_list = [x for x in cell_lines_list if str(x) != 'nan']
                    cell_lines_list.insert(0, 'All')
                    print("Cell lines:", cell_lines_list)
                    cell_lines_widget = widgets.SelectMultiple(options = cell_lines_list, description="Cell line:", value=(cell_lines_list[0],))
                    
                    def update_primary_tissu(tissu):
                        if tissu['new'] != "All":
                            primary_tissu_list = list(set(py_tissues[py_tissues['cell_line'].str.contains('|'.join(tissu['new']), na=False)].primary_disease))
                        else:
                            primary_tissu_list = list(set(py_tissues['primary_disease'].to_list()))
                        print(len(primary_tissu_list))
                        print(primary_tissu_list[0:3])
                        primary_tissu_list.insert(0, 'All')
                        primary_tissu_widget.options = primary_tissu_list
                        primary_tissu_widget.value = (primary_tissu_list[0])
                        
                    def update_cell_lines(cell_lines):
                        if cell_lines['new'] != "All":
                            cell_lines_list_init = list(set(py_tissues[py_tissues['primary_disease'].str.contains('|'.join(cell_lines['new']), na=False)].cell_line_name))
                        else:
                            cell_lines_list_init = list(set(py_tissues['cell_line_name']))
                        cell_lines = [x for x in cell_lines_list_init if str(x) != 'nan']
                        cell_lines.insert(0, 'All')
                        print(len(cell_lines))
                        print(cell_lines[0:3])
                        cell_lines_widget.options = cell_lines
                        cell_lines_widget.value = tuple(cell_lines[0])
                        
                    tissues_widget.observe(update_primary_tissu, 'value')
                    primary_tissu_widget.observe(update_cell_lines, 'value')
                    tissu_selection_button = widgets.Button(description="Run!")
                    
                    display(tissues_widget, primary_tissu_widget, cell_lines_widget, tissu_selection_button)
                    tissu_selection_button.on_click(tissu_selection_button_clicked)

            depmap_query_button = widgets.Button(description="Query!")
            depmap_query_button.on_click(depmap_query_button_clicked)
            display(data_types_widget, depmap_query_button)


        # Output
        global intersection_run
        intersection_run = True
        output = widgets.Output()
        venn_button = widgets.Button(description="Venn diagram", style=dict(button_color='#707070'))
        rra_button = widgets.Button(description="RRA ranking", style=dict(button_color='#707070'))
        genemania_button = widgets.Button(description="Genemania", style=dict(button_color='#707070'))
        enrichr_button = widgets.Button(description="EnrichR", style=dict(button_color='#707070'))
        depmap_button = widgets.Button(description="depmap", style=dict(button_color='#707070'))
        if 'MAGeCK_MLE' in tools_selected:
            params['MAGeCK_MLE']['on'] = True
            mle_fdr = widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05, description='FDR:')
            mle_score=widgets.FloatText(value=0, description='Score cut-off:')
            mle_text=widgets.HTML(value="<b>MAGeCK MLE</b>:")
            mle_order=widgets.ToggleButtons(options=['Greater than score', 'Lower than score'], description='Selection:', name="test")
            display(mle_text)
            mle_box = HBox([mle_fdr, mle_score, mle_order])
            display(mle_box)
            mle_order.observe(mle_order_update, 'value')
            mle_score.observe(mle_score_update, 'value')
            mle_fdr.observe(mle_fdr_update, 'value')
        if 'MAGeCK_RRA' in tools_selected:
            params['MAGeCK_RRA']['on'] = True
            rra_direction = widgets.ToggleButtons(options=['Positive', 'Negative'], description='Direction:')
            rra_fdr = widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05, description='FDR:')
            rra_score=widgets.FloatText(value=0, description='Score cut-off:')
            rra_text=widgets.HTML(value="<b>MAGeCK RRA</b>:")
            rra_order=widgets.ToggleButtons(options=['Greater than score', 'Lower than score'], description='Selection:')
            display(rra_text)
            rra_box = HBox([rra_fdr, rra_direction, rra_score, rra_order])
            display(rra_box)
            rra_direction.observe(rra_direction_update, 'value')
            rra_order.observe(rra_order_update, 'value')
            rra_score.observe(rra_score_update, 'value')
            rra_fdr.observe(rra_fdr_update, 'value')
        if 'BAGEL' in tools_selected:
            params['BAGEL']['on'] = True
            bagel_score=widgets.FloatText(value=0, description='Score cut-off:')
            bagel_text=widgets.HTML(value="<b>BAGEL</b>:")
            bagel_order=widgets.ToggleButtons(options=['Greater than score', 'Lower than score'], description='Selection:')
            display(bagel_text)
            bagel_box = HBox([bagel_score, bagel_order])
            display(bagel_box)
            bagel_order.observe(bagel_order_update, 'value')
            bagel_score.observe(bagel_score_update, 'value')
        if 'CRISPhieRmix' in tools_selected:
            params['CRISPhieRmix']['on'] = True
            CRISPhieRmix_fdr = widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05, description='FDR:')
            CRISPhieRmix_score=widgets.FloatText(value=0, description='Score cut-off:')
            CRISPhieRmix_text=widgets.HTML(value="<b>CRISPhieRmix</b>:")
            CRISPhieRmix_order=widgets.ToggleButtons(options=['Greater than score', 'Lower than score'], description='Selection:')
            display(CRISPhieRmix_text)
            CRISPhieRmix_box = HBox([CRISPhieRmix_fdr, CRISPhieRmix_score, CRISPhieRmix_order])
            display(CRISPhieRmix_box)
            CRISPhieRmix_order.observe(CRISPhieRmix_order_update, 'value')
            CRISPhieRmix_score.observe(CRISPhieRmix_score_update, 'value')
            CRISPhieRmix_fdr.observe(CRISPhieRmix_fdr_update, 'value')
        if 'in_house_method' in tools_selected:
            params['in_house_method']['on'] = True
            in_house_method_direction=widgets.ToggleButtons(options=['Positive', 'Negative'], description='Direction:')
            in_house_method_score=widgets.FloatText(value=0, description='Score cut-off:')
            in_house_method_text=widgets.HTML(value="<b>In-house method</b>:")
            in_house_method_order=widgets.ToggleButtons(options=['Greater than score', 'Lower than score'], description='Selection:')
            display(in_house_method_text)
            in_house_method_box = HBox([in_house_method_direction, in_house_method_score, in_house_method_order])
            display(in_house_method_box)
            in_house_method_order.observe(in_house_method_order_update, 'value')
            in_house_method_score.observe(in_house_method_score_update, 'value')
            in_house_method_direction.observe(in_house_method_direction_update, 'value')
        if 'GSEA-like' in tools_selected:
            params['GSEA_like']['on'] = True
            GSEA_like_fdr = widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05, description='FDR:')
            GSEA_like_score=widgets.FloatText(value=0, description='Score cut-off:')
            GSEA_like_text=widgets.HTML(value="<b>GSEA-like</b>:")
            GSEA_like_order=widgets.ToggleButtons(options=['Greater than score', 'Lower than score'], description='Selection:')
            display(GSEA_like_text)
            GSEA_like_box = HBox([GSEA_like_fdr, GSEA_like_score, GSEA_like_order])
            display(GSEA_like_box)
            GSEA_like_order.observe(GSEA_like_order_update, 'value')
            GSEA_like_score.observe(GSEA_like_score_update, 'value')
            GSEA_like_fdr.observe(GSEA_like_fdr_update, 'value')
        buttons_box = HBox([venn_button, rra_button, genemania_button, enrichr_button, depmap_button])
        display(buttons_box, output)
        venn_button.on_click(venn_button_clicked)
        rra_button.on_click(rra_button_clicked)
        genemania_button.on_click(genemania_button_clicked)
        enrichr_button.on_click(enrichr_button_clicked)
        depmap_button.on_click(depmap_button_clicked)

    _ = interact(run, condition_selected=conditions_widget, tools_selected=tools_widget)



def plot_chart(data, tool, condition, file, x, y, method, column_filter='', cutoff=0, greater=False, element='', column_element='', elem_color='', pass_color='', fail_color='', category_column = '', color_scheme='', reverse=''):
    source = data[tool][condition][file].copy()
    show_plot = True
    if method == 'Cut-off':
        if (source[column_filter].dtype in [np.float64, np.int64]):
            elements = element.split(',')
            if not greater:
                pass_value = "%s < %s" % (column_filter,cutoff)
                fail_value = "%s >= %s" % (column_filter,cutoff)
                source.loc[source[column_filter] < cutoff, 'filter'] = pass_value
                source.loc[source[column_filter] >= cutoff, 'filter'] = fail_value
            else:
                pass_value = "%s >= %s" % (column_filter,cutoff)
                fail_value = "%s < %s" % (column_filter,cutoff)
                source.loc[source[column_filter] >= cutoff, 'filter'] = pass_value
                source.loc[source[column_filter] < cutoff, 'filter'] = fail_value
            source['filter'] = np.where(source[column_element].isin(elements), element, source['filter'])
            domain = [pass_value, fail_value, element]
            range_ = [pass_color, fail_color, elem_color]
        else:
            source['filter'] = 'pass'
            show_plot = False
            print('Choose a column of integer or float type.')
    if show_plot:
        if method == 'Cut-off':
            chart = alt.Chart(source).transform_calculate(
                order="{'%s': 2, '%s': 1, '%s': 0}[datum.filter]" % (element, pass_value, fail_value)
            ).mark_circle().encode(
                x=x,
                y=y,
                tooltip = list(source.columns),
                color=alt.Color('filter', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Pass filter:")),
                order='order:Q'
            ).interactive()
            display(chart)
        else:
            chart = alt.Chart(source).mark_circle().encode(
                x=x,
                y=y,
                tooltip = list(source.columns),
                color=alt.Color(category_column, legend=alt.Legend(title="%s:" % category_column), scale=alt.Scale(scheme=color_scheme, reverse=reverse))
            ).interactive()
            display(chart)


def call_form(tools_available):
    SCHEMES = ["reds", "accent", "redblue", "rainbow"]
    tools = list(tools_available.keys())
    display(widgets.HTML(value="First, choose a <b>tool</b> to browse results:"))
    @interact(
        tool = widgets.Dropdown(options=tools,value=tools[0],description='Tool:',disabled=False))
    def form(tool):
        display(widgets.HTML(value="Choose a <b>condition</b> available for this tool:"))
        conditions = list(tools_available[tool].keys())
        @interact(
            condition = widgets.Dropdown(options=conditions,value=conditions[0],description='Condition:',disabled=False))
        def form(condition):
            display(widgets.HTML(value="Choose a <b>file</b>:"))
            files = list(tools_available[tool][condition].keys())
            @interact(file = widgets.Dropdown(options=files,value=files[0],description='File:',disabled=False))
            def form(file):
                display(widgets.HTML(value="Choose which columns to use as <b>x</b> and <b>y</b> axis:"))
                columns = tools_available[tool][condition][file].columns
                @interact(x = widgets.Dropdown(options=columns,value=columns[0],description='x:',disabled=False),
                          y = widgets.Dropdown(options=columns,value=columns[1],description='y:',disabled=False))
                def form(x , y):
                    display(widgets.HTML(value="Results can be colored using two <b>methods</b>.</br> <ul><li>Cut-off: Choose a <b>column</b> with numerical values on which you want to apply a <b>cut-off</b></li><li>Category: choose a column and a color scheme.</li></ul>"))
                    @interact(method = widgets.ToggleButtons(options=['Cut-off', 'Category'],value='Cut-off', description='Method:'))
                    def form(method):
                        if method == 'Cut-off':
                            display(widgets.HTML(value="Now choose a column with numerical data and a set a <b>cut-off</b> value."))
                            @interact(column_filter = widgets.Dropdown(options=columns,value=columns[0],description='Column:',disabled=False),
                                      cutoff = widgets.FloatText(value=0.0,description='Cut-off:'),
                                      greater = widgets.Checkbox(value=False,description='Greater?',disabled=False,indent=True),
                                      pass_color = widgets.ColorPicker(concise=False,description='Pass color:',value='red',disabled=False),
                                      fail_color = widgets.ColorPicker(concise=False,description='Fail color:',value='gray',disabled=False),
                                      column_element = widgets.Dropdown(options=columns,value=columns[0],description='Column:',disabled=False))
                            def form(column_filter, cutoff, greater, pass_color, fail_color, column_element):
                                display(widgets.HTML(value="It is possible to highlight an <b>element</b> based on a column and value from that column."))
                                first_element = str(tools_available[tool][condition][file][column_element][0])
                                @interact(element = widgets.Text(value=first_element,placeholder='Element:',description='Element:',disabled=False),
                                          elem_color = widgets.ColorPicker(concise=False,description='Element color:',value='blue',disabled=False))
                                def form(element, elem_color):
                                    plot_chart(tools_available, tool, condition, file, x, y, method, column_filter, cutoff, greater, element, column_element, elem_color, pass_color, fail_color)
                        elif method == 'Category':
                            @interact(category_column = widgets.Dropdown(options=columns,value=columns[0],description='Category:',disabled=False),
                                      reverse = widgets.Checkbox(value=False,description='Reverse?',disabled=False,indent=True),
                                      color_scheme = widgets.Dropdown(options=SCHEMES,value=SCHEMES[0],description='Scheme:',disabled=False))
                            def form(category_column, reverse, color_scheme):
                                plot_chart(tools_available, tool, condition, file, x, y, method, category_column=category_column, color_scheme=color_scheme, reverse=reverse)
