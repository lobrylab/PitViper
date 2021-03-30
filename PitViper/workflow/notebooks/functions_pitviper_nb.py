import os
import pandas as pd
import re
import altair as alt
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import requests
import json


class ToolResult:
    TOOLS_PATH = {
        'MAGeCK_MLE': '{tool_path}{directory}/{directory}.gene_summary.txt',
        'MAGeCK_RRA': '{tool_path}{directory}/{directory}.gene_summary.txt',
        'BAGEL': '{tool_path}{directory}/{directory}_BAGEL_output.bf',
        'CRISPhieRmix': '{tool_path}{directory}/{directory}.txt'
        }
    
    SEPARATORS = {
        'MAGeCK_MLE': '\t',
        'MAGeCK_RRA': '\t',
        'BAGEL': '\t',
        'CRISPhieRmix': ','
        }
    def __init__(self, path):
        self.path = path
            
            
    def get_tool_name(self):
        m = re.match("^\.\.\/\.\.\/\w+\/\w+\/(\w+)\/$", self.path)
        if m:
            return m.group(1)
        else:
            print('Failed to get tool name.')
            
            
    def create_comparisons_dict(self):
        comparisons_dict = {}
        self.tool = self.get_tool_name()
        for directory in os.listdir(str(self.path)):
            table_file = self.TOOLS_PATH[self.tool].format(tool_path=self.path, directory=directory)
            comparisons_dict[directory] = {'file': table_file, 
                                           'table': pd.read_csv(table_file, sep=self.SEPARATORS[self.tool])}
        self.comparisons_dict = comparisons_dict

    def init(self):
        self.create_comparisons_dict()
        
        
def create_results_pitviper(results_path):
    results = {}
    for path in results_path:
        res = ToolResult(path=path)
        res.init()
        results[res.tool] = res
    return results

        
def natural_sort(l):
    """Function for natural sorting of list."""
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


def mageck_mle_feature_accros_conditions(mle, baseline, feature, fdr_cutoff, comparisons, treatments):
    
    tables_to_concatenate = []
    for comparison in comparisons:
        table = mle.comparisons_dict[comparison]['table']
        tables_to_concatenate.append(table)
        
    result = pd.concat(tables_to_concatenate, axis=1)
    result = result.iloc[:,~result.columns.duplicated()]
    
    info = result.loc[(result['Gene'] == feature)]
    info = info.filter(regex=r'Gene$|\|beta$|\|z$|\|fdr|\|p-value$')
    info = info.melt(id_vars=('Gene'))
    info[['condition','variable']] = info['variable'].str.split("|", expand=True)
    info = info.pivot(index='condition', columns='variable', values='value')
    info.loc[info['fdr'] < fdr_cutoff, 'significant'] = 'Yes' 
    info.loc[info['fdr'] >= fdr_cutoff, 'significant'] = 'No'
    info.loc[info.index == baseline, 'significant'] = 'Baseline'
    info['feature'] = feature
    info['condition'] = info.index

    domain = ['Yes', 'No', 'Baseline']
    range_ = ['red', 'grey', 'black']

    sort_cols = natural_sort([baseline] + treatments)

    plot = alt.Chart(info).mark_circle(size=60).mark_point(
        filled=True,
        size=100,
        ).encode(
                y='beta',
                x=alt.X('condition:N', sort=sort_cols),
                color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
                tooltip=['feature', 'beta', 'fdr', 'significant'],
        ).properties(
                title=feature + " beta versus baseline (MAGeCK MLE)",
                width=100
        )

    return plot


def mageck_rra_feature_accros_conditions(rra, baseline, feature, fdr_cutoff, comparisons, treatments):
    tables_to_concatenate = []
    for comparison in comparisons:
        table = rra.comparisons_dict[comparison]['table'].copy()
        trt = comparison.split('_vs_')[0]
        new_names = [(i,trt + "_" + i) for i in table.iloc[:, 2:].columns.values]
        table.rename(columns = dict(new_names), inplace=True)
        tables_to_concatenate.append(table)
        
    result = pd.concat(tables_to_concatenate, axis=1)
    result = result.iloc[:,~result.columns.duplicated()]
    
    info = result.loc[(result['id'] == feature)]
    info = info.filter(regex=r'id|.+_neg.+')
    info = info.melt(id_vars=('id'))
    info[['condition','variable']] = info['variable'].str.split("|", expand=True)
    info = info.pivot(index='condition', columns='variable', values='value')
    info.loc[info['fdr'] < fdr_cutoff, 'significant'] = 'Yes' 
    info.loc[info['fdr'] >= fdr_cutoff, 'significant'] = 'No'
    info.loc[info.index == baseline, 'significant'] = 'Baseline'
    info['feature'] = feature
    info['condition'] = info.index
    
    new_row = {'feature':feature, 'condition':baseline, 'significant': 'Baseline', 'lfc':0}
    info = info.append(new_row, ignore_index=True)

    
    domain = ['Yes', 'No', 'Baseline']
    range_ = ['red', 'grey', 'black']

    sort_cols = natural_sort([baseline] + treatments)

    plot = alt.Chart(info).mark_circle(size=60).mark_point(
        filled=True,
        size=100,
        ).encode(
                y='lfc',
                x=alt.X('condition:N', sort=sort_cols),
                color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
                tooltip=['feature', 'lfc', 'fdr', 'significant', 'rank', 'goodsgrna'],
        ).properties(
                title=feature + " lfc versus baseline (MAGeCK RRA)",
                width=100
        )

    
    return plot


def crisphiermix_feature_accros_conditions(crm, baseline, feature, fdr_cutoff, comparisons, treatments):
    tables_to_concatenate = []
    for comparison in comparisons:
        table = crm.comparisons_dict[comparison]['table'].copy()
        trt = comparison.split('_vs_')[0]
        new_names = [(i,trt + "_" + i) for i in table.iloc[:, 1:].columns.values]
        table.rename(columns = dict(new_names), inplace=True)
        tables_to_concatenate.append(table)
        
    result = pd.concat(tables_to_concatenate, axis=1)
    result = result.iloc[:,~result.columns.duplicated()]
         
    info = result.loc[(result['gene'] == feature)]
    info = info.melt(id_vars=('gene'))
    info[['condition','variable']] = info['variable'].str.split("_", expand=True)
    info = info.pivot(index='condition', columns='variable', values='value')
    info.loc[info['locfdr'] < fdr_cutoff, 'significant'] = 'Yes' 
    info.loc[info['locfdr'] >= fdr_cutoff, 'significant'] = 'No'
    info['feature'] = feature
    info['condition'] = info.index
    new_row = {'feature':feature, 'condition':baseline, 'significant': 'Baseline', 'locfdr':None, 'score': 0}
    info = info.append(new_row, ignore_index=True)
    
    domain = ['Yes', 'No', 'Baseline']
    range_ = ['red', 'grey', 'black']

    sort_cols = natural_sort([baseline] + treatments)

    plot = alt.Chart(info).mark_circle(size=60).mark_point(
        filled=True,
        size=100,
        ).encode(
                y='score',
                x=alt.X('condition:N', sort=sort_cols),
                color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
                tooltip=['feature', 'locfdr', 'significant', 'score'],
        ).properties(
                title=feature + " locfc versus baseline (CRISPhieRmix)",
                width=100
        )
    return plot


def bagel_feature_accros_conditions(bag, baseline, feature, fdr_cutoff, comparisons, treatments):
    tables_to_concatenate = []
    for comparison in comparisons:
        table = bag.comparisons_dict[comparison]['table'].copy()
        trt = comparison.split('_vs_')[0]
        new_names = [(i,trt + "_" + i) for i in table.iloc[:, 1:].columns.values]
        table.rename(columns = dict(new_names), inplace=True)
        tables_to_concatenate.append(table)
        
    result = pd.concat(tables_to_concatenate, axis=1)
    result = result.iloc[:,~result.columns.duplicated()]
         
    info = result.loc[(result['GENE'] == feature)]
    info = info.melt(id_vars=('GENE'))
    info[['condition','variable']] = info['variable'].str.split("_", expand=True)
    info = info.pivot(index='condition', columns='variable', values='value')
    info.loc[info['BF'] < 0, 'significant'] = 'No' 
    info.loc[info['BF'] >= 0, 'significant'] = 'Yes'
    info['feature'] = feature
    info['condition'] = info.index
    new_row = {'feature':feature, 'condition':baseline, 'significant': 'Baseline', 'BF':0, 'STD': None, 'NumObs': None}
    info = info.append(new_row, ignore_index=True)
    
    domain = ['Yes', 'No', 'Baseline']
    range_ = ['red', 'grey', 'black']

    sort_cols = natural_sort([baseline] + treatments)

    plot = alt.Chart(info).mark_circle(size=60).mark_point(
        filled=True,
        size=100,
        ).encode(
                y='BF',
                x=alt.X('condition:N', sort=sort_cols),
                color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
                tooltip=['feature', 'BF', 'significant', 'STD', 'NumObs'],
        ).properties(
                title=feature + " BF versus baseline (BAGEL)",
                width=100
        )
    return plot
    
    
def feature_accros_conditions(pitviper_res):
    conditions = [{'baseline':condition.split('_vs_')[1] , 'treatment':condition.split('_vs_')[0]} for condition in pitviper_res.comparisons_dict.keys()]
    @interact(baseline=set([condition['baseline'] for condition in conditions]), feature=widgets.Text(value='MYC', placeholder='Feature to show...', description='Feature:'))
    def feature_accros_conditions(baseline, feature):
        @interact(fdr_cutoff=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05))
        def feature_accros_conditions(fdr_cutoff):
            treatments = [condition['treatment'] for condition in conditions if condition['baseline'] == baseline]
            comparisons = [condition['treatment']+'_vs_'+baseline for condition in conditions if condition['baseline'] == baseline]
            if pitviper_res.tool == 'MAGeCK_MLE':
                plot = mageck_mle_feature_accros_conditions(pitviper_res, baseline, feature, fdr_cutoff, comparisons, treatments)
            if pitviper_res.tool == 'MAGeCK_RRA':
                plot = mageck_rra_feature_accros_conditions(pitviper_res, baseline, feature, fdr_cutoff, comparisons, treatments)
            if pitviper_res.tool == 'CRISPhieRmix':
                plot = crisphiermix_feature_accros_conditions(pitviper_res, baseline, feature, fdr_cutoff, comparisons, treatments)
            if pitviper_res.tool == 'BAGEL':
                plot = bagel_feature_accros_conditions(pitviper_res, baseline, feature, fdr_cutoff, comparisons, treatments)   
            return plot
        
        
        
        
        
        
        
        
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


def enrichmentBarPlot(source, n, description):
    source = source.sort_values(by=['Combined score'], ascending=False).head(n)

    domain = [source['Adjusted p-value'].min(), source['Adjusted p-value'].max()]
    range_ = ['darkred', 'white']

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
        title=description,
    )


    chart = (bars).properties(height=15*n, width=500)
    return chart


def enrichmentCirclePlot(source, n):
    source = source.sort_values(by=['Combined score'], ascending=False).head(n)

    domain = [source['Adjusted p-value'].min(), source['Adjusted p-value'].max()]
    range_ = ['darkred', 'white']
    
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
                 'Old adjusted p-value',],
    ).interactive()
    
    chart = (chart).properties(height=300, width=900)
    
    return chart