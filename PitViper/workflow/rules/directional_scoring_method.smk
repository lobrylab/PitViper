



rule counts_filtering:
    input:
        normalized=config['normalized_count_table'],
        raw=config['count_table_file']
    output:
        normalized_filtered_counts = "results/" + config['token'] + "/normalized.filtered.counts.txt",
        raw_filtered_counts = "results/" + config['token'] + "/raw.filtered.counts.txt"
    params:
        design=config['tsv_file'],
        threshold=config['counts_threshold']
    log:
        "logs/" + config['token'] + "/counts_filtering/" + "counts_filtering.log"
    script:
        "../../workflow/scripts/counts_filtering.py"
