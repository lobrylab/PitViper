rule counts_filtering:
    """ Filter counts based on a threshold. """
    input:
        normalized=config['normalized_count_table'],
        raw=config['count_table_file']
    output:
        normalized_filtered_counts = f"results/{config['token']}/normalized.filtered.counts.txt",
        raw_filtered_counts = f"results/{config['token']}/raw.filtered.counts.txt"
    params:
        design=config['tsv_file'],
        threshold=config['counts_threshold']
    log:
        f"logs/{config['token']}/counts_filtering/counts_filtering.log"
    message:
        "Filtering counts with threshold {params.threshold} from {input.raw} and {input.normalized}."
    script:
        "../../workflow/scripts/counts_filtering.py"
