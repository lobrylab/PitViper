rule MAGeCK_counts_normalize:
    """Normalizes the counts table using MAGeCK count."""
    input:
        config['count_table_file']
    output:
        config['normalized_count_table']
    params:
        name = f"resources/{config['token']}/screen",
    log:
        f"logs/{config['token']}/MAGeCK_counts_normalize.log"
    message:
        "Normalizing counts table: {input} to {output}."
    shell:
        "mageck count \
            -k {input} \
            -n {params.name} \
            --norm-method total > {log} 2>&1"