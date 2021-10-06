import os
import pandas as pd
import re
import altair as alt
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import requests
import json
import natsort as ns
import numpy as np
import yaml
from mc4.algorithm import mc4_aggregator
import os.path
from os import path


def open_yaml(yml):
    """Open a YAML file and return it's content."""
    with open(yml, "r") as stream:
        try:
            content = yaml.safe_load(stream)
            return content
        except yaml.YAMLError as exc:
            print(exc)
            

def setup_step_1(token):
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
            
    return (results_directory, tools_available)
            
            
            
            
def show_mapping_qc(token):
    path_qc = "./results/%s/screen.countsummary.txt" % token
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
    path_qc = content['count_table_file']
    if not path.exists(path_qc):
        print("No count file to show.")
        return 0
    table = pd.read_csv(path_qc, sep='\t')
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
        new_row = {'pathway':gene, 'condition':control, 'significant': 'Baseline', 'pval':1, 'padj':1 ,'ES': 0, 'NES': 0, 'nMoreExtreme':None, 'size':None}
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
                    data = data.rename(columns=lambda x: re.sub('\w+\|','',x))
                    tables_list.append(data)
                else:
                    tables_list.append(data)
    result = pd.concat(tables_list)
    return result


def MAGeCK_MLE_results(results_directory, tools_available):
    tool = "MAGeCK_MLE"
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
    
    def on_button_clicked(b):
        chart = alt.Chart(source).mark_circle(size=60).encode(
            x=alt.X('default_rank:Q', axis=alt.Axis(title='Rank')),
            y=alt.Y(treatment + '|beta:Q'),
            tooltip=['Gene', 'sgRNA', treatment + '|beta', treatment + '|fdr', 'significant'],
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
    
    
    
    
def MAGeCK_RRA_snake_plot(comparison, fdr_cutoff, non_sig, sig, results_directory, tools_available):
    tool = "MAGeCK_RRA"
    print(tool)
    print(comparison)
    print(fdr_cutoff)
    treatment, control = comparison.split("_vs_")
    source = MAGeCK_RRA_data(comparison = comparison, control = "", tool = tool, results_directory=results_directory, tools_available=tools_available)

    source['default_rank'] = source['neg|lfc'].rank()
    source.loc[source['neg|fdr'] < 0.05, 'significant'] = 'Yes' 
    source.loc[source['neg|fdr'] >= 0.05, 'significant'] = 'No'
    domain = ['Yes', 'No']
    range_ = [sig, non_sig]
    
    def on_button_clicked(b):
        chart = alt.Chart(source).mark_circle(size=60).encode(
            x=alt.X('default_rank:Q', axis=alt.Axis(title='Rank')),
            y=alt.Y('neg|lfc:Q'),
            tooltip=['id', 'num', 'neg|lfc', 'neg|fdr', 'significant', 'neg|rank'],
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
    source.loc[source['FDR'] < 0.05, 'significant'] = 'Yes' 
    source.loc[source['FDR'] >= 0.05, 'significant'] = 'No'
    domain = ['Yes', 'No']
    range_ = [sig, non_sig]
    
    def on_button_clicked(b):
        chart = alt.Chart(source).mark_circle(size=60).encode(
            x=alt.X('default_rank:Q', axis=alt.Axis(title='Rank')),
            y=alt.Y('locfdr:Q'),
            tooltip=['gene', 'locfdr', 'score', 'FDR', 'significant'],
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
    source.loc[abs(source['score']) > 0.5, 'significant'] = 'Yes' 
    source.loc[abs(source['score']) <= 0.5, 'significant'] = 'No'
    domain = ['Yes', 'No']
    range_ = [sig, non_sig]
    
    def on_button_clicked(b):
        chart = alt.Chart(source).mark_circle(size=60).encode(
            x=alt.X('default_rank:Q', axis=alt.Axis(title='Rank')),
            y=alt.Y('score:Q'),
            tooltip=['Gene', 'up', 'down', 'n', 'significant', 'score'],
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
    
    
    
    
def GSEA_like_snake_plot(comparison, fdr_cutoff, non_sig, sig):
    tool = "GSEA-like"
    print(tool)
    print(comparison)
    print(fdr_cutoff)
    treatment, control = comparison.split("_vs_")
    source = GSEA_like_data(comparison = comparison, control = "", tool = tool, results_directory=results_directory, tools_available=tools_available)
    source['default_rank'] = source['NES'].rank(method='dense')
    source.loc[abs(source['pval']) < fdr_cutoff, 'significant'] = 'Yes' 
    source.loc[abs(source['pval']) >= fdr_cutoff, 'significant'] = 'No'
    domain = ['Yes', 'No']
    range_ = [sig, non_sig]
    
    def on_button_clicked(b):
        chart = alt.Chart(source).mark_circle(size=60).encode(
            x=alt.X('default_rank:Q', axis=alt.Axis(title='Rank')),
            y=alt.Y('NES:Q'),
            tooltip=['pathway', 'pval', 'padj', 'ES', 'NES','significant', 'nMoreExtreme', 'size'],
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
            tooltip=['GENE', 'BF', 'STD', 'significant'],
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
    @interact(tools=widgets.Dropdown(options=set(["CRISPhieRmix", "GSEA-like", "in_house_method", "MAGeCK_MLE", "MAGeCK_RRA", "BAGEL"]), description='Tool:', disabled=False))
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