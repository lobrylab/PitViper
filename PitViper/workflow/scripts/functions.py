import re
import os
import pandas as pd
import time
import altair as alt
from functools import reduce

from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets



def natural_sort(l):
    """Function for natural sorting of list."""
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)



def mageck_mle_viz(mageck_mle_outputs):
    def mageck_mle_inputs(mageck_mle_outputs):
        mageck_mle = {}
        for directory in os.listdir(mageck_mle_outputs):
            table_file = mageck_mle_outputs+directory+'/'+directory+'.gene_summary.txt'
            mageck_mle[directory] = {'file': table_file, 'table': pd.read_csv(table_file, sep="\t")}
        return mageck_mle
    

    def get_gene_info(table, gene, col):
        gene_info = table.loc[(table[col] == gene)]
        return gene_info


    def get_gene_accross_conditions_plot(genes_summary, conditions, gene, fdr_cutoff, baseline):
        rows = []
        for condition in conditions:
            info = get_gene_info(genes_summary, gene, 'Gene')
            beta = '{condition}|beta'.format(condition=condition)
            fdr = '{condition}|fdr'.format(condition=condition)
            info = info[["Gene", fdr, beta]]
            info['condition'] = condition
            info.columns = ['Gene', 'fdr', 'beta', 'condition']
            rows.append(info)


        result = pd.concat(rows)
        result.loc[result['fdr'] < fdr_cutoff, 'significant'] = 'Yes' 
        result.loc[result['fdr'] >= fdr_cutoff, 'significant'] = 'No'
        result.loc[result['condition'] == baseline, 'significant'] = 'Baseline'

        chart = alt.Chart(result)

        sort_cols = natural_sort(conditions)

        domain = ['Yes', 'No', 'Baseline']
        range_ = ['red', 'grey', 'black']

        plot = alt.Chart(result).mark_circle(size=60).mark_point(
            filled=True,
            size=100,
            ).encode(
                    y='beta',
                    x=alt.X('condition:N', sort=sort_cols),
                    color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
                    tooltip=['Gene', 'beta', 'fdr', 'significant'],
            ).properties(
                    title=gene + " beta versus baseline (MAGeCK MLE)",
                    width=100
            )
        return plot
    
    
    def columns_to_drop_mageck_mle(table, baseline):
        z = '{}|z'.format(baseline)
        p_value = '{}|p-value'.format(baseline)
        waldp = '{}|wald-p-value'.format(baseline)
        waldfdr = '{}|wald-fdr'.format(baseline)
        cols = [z, p_value, waldp, waldfdr]
        return cols
    
    
    def displayTable(table, conditions):
        print("\nTable's fdr cutoff for given condition:")
        @interact(fdr=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05),
                  conditions=conditions)
        def displayTable(fdr, conditions):
            return table[table['{conditions}|fdr'.format(conditions=conditions)] < fdr]


    def mageck_mle_results(mageck_mle, control):
        to_concat = []
        conditions_plot = []
        conditions_plot.append(control)
        for condition in mageck_mle.keys():
            trt, con = condition.split("_vs_")
            if con == control:
                conditions_plot.append(trt)
                table = mageck_mle[condition]['table']
                table = table.drop(columns_to_drop_mageck_mle(table, control), axis=1)
                to_concat.append(table)

        result = pd.concat(to_concat, axis=1)
        result = result.iloc[:,~result.columns.duplicated()]
        
        return (result, conditions_plot)
    
    
    mageck_mle = mageck_mle_inputs(mageck_mle_outputs)
    print("Chart parameters :")
    @interact(feature=list(mageck_mle[list(mageck_mle.keys())[0]]['table']['Gene'].values)[0], 
              fdr=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05),
              control=list(set([i.split("_vs_")[1] for i in mageck_mle.keys()])))
    def mageck_mle_widgets(feature, fdr, control):
        time.sleep(1)
        result = mageck_mle_results(mageck_mle, control)
        table = result[0]
        conditions = result[1]
        displayTable(table, conditions)
        chart = get_gene_accross_conditions_plot(table, conditions, feature, fdr, control)
        return chart
    
    
    
    
    
def mageck_rra_viz(mageck_rra_outputs):
    def mageck_rra_inputs(mageck_rra_outputs):
        mageck_rra = {}
        for directory in os.listdir(mageck_rra_outputs):
            table_file = mageck_rra_outputs+directory+'/'+directory+'.gene_summary.txt'
            mageck_rra[directory] = {'file': table_file, 'table': pd.read_csv(table_file, sep="\t")}
        return mageck_rra

    def displayTable(table, conditions):
        print("\nTable's fdr cutoff for given condition:")
        @interact(fdr=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05),
                  conditions=conditions)
        def displayTable(fdr, conditions):
            return table[table['{conditions}_neg|fdr'.format(conditions=conditions)] < fdr]

    def mageck_rra_results(mageck_rra, control):
        tables = []
        conditions_plot = []
        for condition in mageck_rra.keys():
            trt, con = condition.split("_vs_")
            if con == control:
                conditions_plot.append(trt)
                table = mageck_rra_inputs(mageck_rra_outputs)[condition]['table']
                new_names = [(i,trt + "_" + i) for i in table.iloc[:, 2:].columns.values]
                table.rename(columns = dict(new_names), inplace=True)
                tables.append(table)

        big_table = reduce(lambda  left,right: pd.merge(left,right,on=['id'],
                                                    how='outer'), tables)

        return (big_table, conditions_plot)
    
    def mageck_rra_plot(genes_summary, conditions, gene, fdr_cutoff, baseline):
        rows = []
        for condition in conditions:
            if condition != baseline:
                info = genes_summary.loc[(genes_summary['id'] == gene)]
                lfc = '{trt}_neg|lfc'.format(trt=condition, con=baseline)
                fdr = '{trt}_neg|fdr'.format(trt=condition, con=baseline)
                info = info[["id", fdr, lfc]]
                info['condition'] = condition
                info.columns = ['id', 'fdr', 'lfc', 'condition']
                rows.append(info)

        result = pd.concat(rows)
        result.loc[result['fdr'] < fdr_cutoff, 'significant'] = 'True' 
        result.loc[result['fdr'] >= fdr_cutoff, 'significant'] = 'False'
        result.loc[result['condition'] == baseline, 'significant'] = 'Baseline'

        chart = alt.Chart(result)

        sort_cols = natural_sort([cond.split("_vs_")[0] for cond in conditions])

        domain = ['True', 'False', 'Baseline']
        range_ = ['red', 'grey', 'black']

        plot = alt.Chart(result).mark_circle(size=60).mark_point(
            filled=True,
            size=100,
            ).encode(
                    y='lfc',
                    x=alt.X('condition:N', sort=sort_cols),
                    color=alt.Color('significant', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Significativity:")),
                    tooltip=['id', 'lfc', 'fdr', 'significant'],
            ).properties(
                    title=gene + " beta versus baseline (MAGeCK RRA)",
                    width=100
            )

        return plot

    mageck_rra = mageck_rra_inputs(mageck_rra_outputs)
    print("Chart parameters :")
    @interact(feature=list(mageck_rra[list(mageck_rra.keys())[0]]['table']['id'].values)[0], 
              fdr=widgets.FloatSlider(min=0.0, max=1.0, step=0.01, value=0.05),
              control=list(set([i.split("_vs_")[1] for i in mageck_rra.keys()])))
    def mageck_rra_widgets(feature, fdr, control):
        time.sleep(1)
        result = mageck_rra_results(mageck_rra, control)
        table = result[0]
        conditions = result[1]
        displayTable(table, conditions)
        chart = mageck_rra_plot(table, conditions, feature, fdr, control)
        return chart