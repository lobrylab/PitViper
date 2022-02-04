



rule counts_filtering:
    input:
        normalized=config['normalized_count_table'],
        raw=config['count_table_file']
    output:
        normalized_filtered_counts = "results/" + config['token'] + "/normalized.filtered.counts.txt",
        raw_filtered_counts = "results/" + config['token'] + "/raw.filtered.counts.txt"
    params:
        config['tsv_file'],
    script:
        "../../workflow/scripts/counts_filtering.py"
