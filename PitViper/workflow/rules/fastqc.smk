

def getFastqcInput(wildcards):
    sample_sheet = pd.read_csv(config['tsv_file'], sep="\t")
    # If column 'fastq' is present, use it
    path_to_files = []
    if 'fastq' in sample_sheet.columns:
        path_to_files = sample_sheet['fastq'].tolist()
    # If column 'bam' is present, use it
    elif 'bam' in sample_sheet.columns:
        path_to_files = sample_sheet['bam'].tolist()
    return path_to_files


def getFastqcInputParams(wildcards):
    sample_sheet = pd.read_csv(config['tsv_file'], sep="\t")
    # If column 'fastq' is present, use it
    path_to_files = []
    if 'fastq' in sample_sheet.columns:
        path_to_files = sample_sheet['fastq'].tolist()
    # If column 'bam' is present, use it
    elif 'bam' in sample_sheet.columns:
        path_to_files = sample_sheet['bam'].tolist()
    # Create a string with files separated by space
    return " ".join(path_to_files)
    

rule fastqc:
    input:
        getFastqcInput
    output:
        directory(f"results/{config['token']}/fastqc")
    params:
        input=getFastqcInputParams
    log:
        f"logs/{config['token']}/fastqc.log"
    shell:
        """
        mkdir -p {output} && \
        fastqc -o {output} {params.input} &> {log}
        """