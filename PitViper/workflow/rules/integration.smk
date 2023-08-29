rule genes_integration:
    """ Integrates the results of the different tools to generate a report. """
    input:
        generatedResults
    output:
        f"results/{config['token']}/Report.ipynb"
    params:
        config['token']
    log:
        notebook=f"results/{config['token']}/Report.ipynb"
    notebook:
        "../notebooks/Report_template.py.ipynb"
