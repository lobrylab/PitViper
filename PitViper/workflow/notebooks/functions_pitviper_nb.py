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
import os.path
from os import path

alt.data_transformers.disable_max_rows()

class ToolResult:
    TOOLS_RESULT_PATH = {
        'MAGeCK_MLE': '{tool_path}{directory}/{directory}.gene_summary.txt',
        'MAGeCK_RRA': '{tool_path}{directory}/{directory}.gene_summary.txt',
        'BAGEL': '{tool_path}{directory}/{directory}_BAGEL_output.bf',
        'CRISPhieRmix': '{tool_path}{directory}/{directory}.txt'
        }
    
    TOOLS_sgRNA_PATH = {
        'MAGeCK_MLE': '{tool_path}{directory}/{directory}.sgrna_summary.txt',
        'MAGeCK_RRA': '{tool_path}{directory}/{directory}.sgrna_summary.txt'}
    
    DESEQ2_PATH = '{tool_path}{directory}/{directory}_DESeq2_results.txt'
    
    SEPARATORS = {
        'MAGeCK_MLE': '\t',
        'MAGeCK_RRA': '\t',
        'BAGEL': '\t',
        'CRISPhieRmix': ','
        }
    def __init__(self, path):
        self.path = path
            
            
    def get_tool_name(self):
        m = re.match("\w+\/\w+\/(\w+)\/$", self.path)
        if m:
            return m.group(1)
        else:
            print('Failed to get tool name.')
            
            
    def create_comparisons_dict(self):
        comparisons_dict = {}
        self.tool = self.get_tool_name()
        for directory in os.listdir(str(self.path)):
            table_file = self.TOOLS_RESULT_PATH[self.tool].format(tool_path=self.path, directory=directory)
            comparisons_dict[directory] = {'file': table_file, 
                                           'table': pd.read_csv(table_file, sep=self.SEPARATORS[self.tool])}
            if self.tool in ['MAGeCK_MLE', 'MAGeCK_RRA']:
                sgRNA_table = self.TOOLS_sgRNA_PATH[self.tool].format(tool_path=self.path, directory=directory)
                comparisons_dict[directory]['sgRNA'] = pd.read_csv(sgRNA_table, sep=self.SEPARATORS[self.tool])
#             if self.tool in ['CRISPhieRmix']:
#                 DESeq2_results = self.DESEQ2_PATH.format(tool_path=self.path, directory=directory)
#                 comparisons_dict[directory]['DESeq2'] = pd.read_csv(DESeq2_results, sep=self.SEPARATORS[self.tool])
                    
        self.comparisons_dict = comparisons_dict

    def init(self):
        self.create_comparisons_dict()
        
        
        
        
def show_read_count_distribution(token, width=800, height=400):
    path_qc = "./results/{token}/screen.count.txt".format(token=token)
    if not path.exists(path_qc):
        print("No count fileto show.")
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
    ).properties(width=width, height=height)#.interactive()
    
    return chart


        
def show_mapping_qc(token):
    path_qc = "./results/{token}/screen.countsummary.txt".format(token=token)
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




def gene_selection_table(results):
    TOOLS = [str(tool) for tool in results.keys()]
    @interact(tool=TOOLS)
    def gene_selection_table(tool):
        tool_res = results[tool]
        conditions = [{'baseline':condition.split('_vs_')[1] , 'treatment':condition.split('_vs_')[0]} for condition in tool_res.comparisons_dict.keys()]
        @interact(baseline=set([condition['baseline'] for condition in conditions]))
        def gene_selection_table(baseline):
            treatments = [condition['treatment'] for condition in conditions if condition['baseline'] == baseline]
            @interact(treatment=treatments)
            def gene_selection_table(treatment):
                print('Tool:', tool)
                print('Baseline:', baseline)
                print('Treatment', treatment)

                def on_button_clicked(b):

                    info = tool_res.comparisons_dict[treatment+'_vs_'+baseline]['table']
                    info = info.drop(info.filter(like=baseline).columns,axis=1)
                    
                    with output:
                        display(info)

                button = widgets.Button(description="Show datatable")
                output = widgets.Output()

                display(button, output)

                button.on_click(on_button_clicked)



def sgRNA_accross_conditions(result, width=500, height=250):
    conditions = [{'baseline':condition.split('_vs_')[1] , 'treatment':condition.split('_vs_')[0]} for condition in result.comparisons_dict.keys()]
    baselines=set([condition['baseline'] for condition in conditions])
    features=widgets.Text(value='', placeholder='Feature to show...', description='Feature:')
    @interact(feature=features, baseline=baselines)
    def sgRNA_accross_conditions(feature, baseline):
        treatments = [condition['treatment'] for condition in conditions if condition['baseline'] == baseline]
        @interact(treatment=treatments)
        def sgRNA_accross_conditions(treatment):
            print('Feature:', feature)
            print('Baseline', baseline)
            print('Treatment', treatment)
            if feature != '':
                tables_to_concatenate = []
                comp = treatment + '_vs_' + baseline
                table = result.comparisons_dict[comp]['sgRNA'].copy()
                new_names = [(i,treatment + "|" + i) for i in table.iloc[:, 2:].columns.values]
                table.rename(columns = dict(new_names), inplace=True)
                tables_to_concatenate.append(table)

                def on_button_clicked(b):
                    info = pd.concat(tables_to_concatenate, axis=1)
                    info = info.iloc[:,~info.columns.duplicated()]

                    info = info.loc[info['Gene'] == feature]

                    info = info.filter(regex=r'sgrna$|\|control_count$|\|treatment_count$')

                    info[[baseline+'-1',baseline+'-2', baseline+'-3']] = info[treatment + '|control_count'].str.split("/",expand=True,)

                    info[[treatment+'-1',treatment+'-2', treatment+'-3']] = info[treatment + '|treatment_count'].str.split("/",expand=True,)

                    info = info.drop(columns=[treatment+'|control_count', treatment+'|treatment_count'])

                    info = info.melt(id_vars=('sgrna'))


                    sort_cols = natural_sort([baseline] + treatments)


                    chart = alt.Chart(info).mark_line().encode(
                        x=alt.X('variable', axis=alt.Axis(title='Condition'), sort=sort_cols),
                        y=alt.X('value:Q', axis=alt.Axis(title='Count'), sort="ascending"),
                        color='sgrna',
                        tooltip=['sgrna', 'value']
                    ).properties(width=width, height=height, title=feature + " sgRNAs read count by condition")
                    with output:
                        display(chart)

                button = widgets.Button(description="Show sgRNA counts")
                output = widgets.Output()

                display(button, output)

                button.on_click(on_button_clicked)
            else: return "No feature to show."


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
                title=feature + " BF versus baseline" (
    ),
                width=100
        )
    return plot
    
    
def feature_accros_conditions(pitviper_res):
    conditions = [{'baseline':condition.split('_vs_')[1] , 'treatment':condition.split('_vs_')[0]} for condition in pitviper_res.comparisons_dict.keys()]
    @interact(baseline=set([condition['baseline'] for condition in conditions]), feature=widgets.Text(value='', placeholder='Feature to show...', description='Feature:'))
    def feature_accros_conditions(baseline, feature):
        @interact(fdr_cutoff=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05))
        def feature_accros_conditions(fdr_cutoff):
            if feature != '':
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
            else: return "no feature to show."
        
        
        
        
        
        
        
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


def enrichmentBarPlot(source, n, description, col_1, col_2):
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
        title=description,
    )


    chart = (bars).properties(height=15*n, width=500)
    return chart


def enrichmentCirclePlot(source, n, description, col_1, col_2):
    if n == 'max':
        n = int(len(source.index))
        
    source = source.sort_values(by=['Combined score'], ascending=False).head(n)

    source['n_overlap'] = source['Overlapping genes'].str.len()
        
    domain = [0, 0.05, 1]
    range_ = [col_1, "grey", "grey"]
    
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
                 'n_overlap'],
        size='n_overlap',
    ).properties(
        title=description
    ).configure_axisY(
        titleAngle=0, 
        titleY=-10,
        titleX=-60,
        labelLimit=1000
    ).interactive()
    
    chart = (chart).properties(height=20*n, width=500)
    
    return chart






        
def enrichr_plots(pitviper_res):
    BASES = open("workflow/notebooks/enrichr_list.txt", "r").readlines()
    TOOLS = [tool for tool in pitviper_res.keys()]
    @interact(tool=TOOLS)
    def enrichr_plots(tool):
        tool_res = pitviper_res[tool]
        conditions = [{'baseline':condition.split('_vs_')[1] , 'treatment':condition.split('_vs_')[0]} for condition in tool_res.comparisons_dict.keys()]
        @interact(description=widgets.Text(value='My gene list', placeholder='Description', description='Description:'), base=BASES, baseline=set([condition['baseline'] for condition in conditions]))
        def enrichr_plots(description, base, baseline):
            treatments = [condition['treatment'] for condition in conditions if condition['baseline'] == baseline]
            @interact(col_2=widgets.ColorPicker(concise=False, description='Top color', value='blue', disabled=False), col_1=widgets.ColorPicker(concise=False, description='Bottom color', value='red', disabled=False), plot_type=['Circle', 'Bar'], size=[5, 10, 20, 50, 100, 200, 'max'], treatment=treatments, fdr_cutoff=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05), score_cutoff=widgets.Text(value='0', placeholder='0', description='Score cut-off:'))
            def enrichr_plots(treatment, score_cutoff, fdr_cutoff, size, plot_type, col_2, col_1):
                comparisons = [condition['treatment']+'_vs_'+baseline for condition in conditions if condition['baseline'] == baseline]
                print('Description:', description)
                print('Baseline:', baseline)
                print('Treatment:', treatment)
                print('FDR cut-off:', fdr_cutoff)
                print('Score cut-off', score_cutoff)
                print('Tool:', tool)
                print('Gene set library:', base)

                score_cutoff = float(score_cutoff)
                
                def on_button_clicked(b):
                    

                    if tool == 'MAGeCK_MLE':
                        info = tool_res.comparisons_dict[treatment+'_vs_'+baseline]['table']
                        info = info.loc[info[treatment+'|fdr'] < fdr_cutoff]
                        if score_cutoff < 0:
                            info = info.loc[info[treatment+'|beta'] < score_cutoff]
                        elif score_cutoff > 0:
                            info = info.loc[info[treatment+'|beta'] > score_cutoff]
                        genes = info['Gene']

                    if tool == 'MAGeCK_RRA':
                        info = tool_res.comparisons_dict[treatment+'_vs_'+baseline]['table']
                        print(info.columns)
                        info = info.loc[info['neg|fdr'] < fdr_cutoff]
                        genes = info['id']
                        
                    print("Size (gene set):", len(genes))

                    enrichr_res = getEnrichrResults(genes, description, base)
                    table = createEnrichrTable(enrichr_res)
                    if plot_type == 'Bar':
                        chart = enrichmentBarPlot(table, size, description, col_1, col_2)
                    else:
                        chart = enrichmentCirclePlot(table, size, description, col_1, col_2)
                    with output:
                        display(chart)


                button = widgets.Button(description="Show EnrichR results")
                output = widgets.Output()

                display(button, output)

                button.on_click(on_button_clicked)
                
                
                
                
def test_plotting(data, width=800, height=400):
    conditions = [{'baseline':condition.split('_vs_')[1] , 'treatment':condition.split('_vs_')[0]} for condition in data.comparisons_dict.keys()]
    baselines=set([condition['baseline'] for condition in conditions])
    @interact(baseline=baselines)
    def test_plotting(baseline):
        treatments = [condition['treatment'] for condition in conditions if condition['baseline'] == baseline]
        non_sig=widgets.ColorPicker(concise=False, description='Non-significant', value='#949494', disabled=False)
        sig=widgets.ColorPicker(concise=False, description='Significant', value='red', disabled=False)
        @interact(treatment=treatments, non_sig=non_sig, sig=sig)
        def test_plotting(treatment, non_sig, sig):
            print('Baseline:', baseline)
            print('Treatment:', treatment)
            comp = treatment + '_vs_' + baseline
            source = data.comparisons_dict[comp]['table']

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
                    color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:"))
                ).properties(width=width, height=height).interactive()

                line = alt.Chart(pd.DataFrame({'y': [0]})).mark_rule().encode(y='y')

                chart = (chart + line)
                with output:
                    display(chart)
            
            button = widgets.Button(description="Show plot")
            output = widgets.Output()

            display(button, output)

            button.on_click(on_button_clicked)
          
