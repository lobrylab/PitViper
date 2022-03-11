<img src="PitViper/docs/logo/pitviper_remasteredv2.png" alt="alt text" width="500" height="175">

> **A pipeline to analyze functional screening data from shRNA, CRISPR/Cas9 and CRISPR/dCas9 experiments.**

## Introduction

PitViper is intended to facilitate analysis of functional screening data from various experiments (shRNA, CRISPR/Cas9 or CRISPR/dCas9).

The pipeline is built with [`Snakemake`](https://snakemake.readthedocs.io/en/stable/), a workflow management system to create reproducible and scalable data analysis, [`Flask`](https://flask.palletsprojects.com/en/2.0.x/), a lightweight web framework and [`Jupyter`](https://jupyter.org/), a web application for creating and sharing computational documents.


##### Table of contents

- [Prerequisites](#prerequisites)

- [Installation](#installation)

- [Inputs](#inputs)

  - [Starting from raw FASTQ or BAM files](#starting-from-raw-fastq-files)

  - [Starting from count matrix](#starting-from-count-matrix)


### Prerequisites

To install and run PitViper, [`Conda`](https://docs.conda.io/en/latest/) and [`Git`](https://git-scm.com/) need to be available from commande-line.

Once Conda is installed and in order to speed-up installation process, [`Mamba`](https://github.com/mamba-org/mamba) is needed:

```bash
$ conda install -c conda-forge mamba
```

Then, clone PitViper repository using `Git`:

```bash
$ git clone https://github.com/PaulArthurM/PitViper.git  # Clone PitViper reposity in ~/PitViper/

$ cd PitViper/PitViper  # Change working directory to PitViper root directory
```

### Installation

Now, install PitViper dependancies. To facilitate dependancy management, a Conda YAML file containing all dependancies have been created and can be used along with `Mamba` to automatically install all dependencies: `pitviper_env.yaml`.

Furthermore, `install_PitViper_env.sh` is a bash script created to perform this step in one command:

```bash
$ ./install_PitViper_env.sh
```

`./install_PitViper_env.sh` first create a Conda environment call `pitviper_env` that should be visible once installation is over, as follow:

```bash
$ conda env list
  # conda environments:
  #
  base                  *  /home/paularthur/miniconda3
  pitviper_env             /home/paularthur/miniconda3/envs/pitviper_env
```

Once `pitviper_env` is created, you can run `run.sh` script from command-line, it should run a server with `Flask` in background and open PitViper GUI in your default web browser:

```bash
$ ./run.sh
```


<img src="PitViper/docs/PitViper.png" alt="alt text">


### Inputs

Prior to use PitViper, we need to introduce it's input files.

PitViper allow to start an analysis from raw data files such as FASTQ, already aligned BAM file or a count matrix. Depending of that, input files will be different.

#### Starting from raw FASTQ files

You will need:

1. Path to FATSQ files on system
2. A library file with three comma-separated columns without header: shRNA ID, shRNA sequence, targer element.
3. A design matrix that summary your conditions, their replicates and associated FASTQ files.

##### Design matrix

Let say that you have two conditions, A (control) and B (treatment), with 3 replicates for each. The associated matrix should look as below:

| condition | replicate | fastq                 | order |
|-----------|-----------|-----------------------|-------|
| A         | A_1       | /path/to/A_rep1.fastq | 0     |
| A         | A_2       | /path/to/A_rep2.fastq | 0     |
| A         | A_3       | /path/to/A_rep3.fastq | 0     |
| B         | B_1       | /path/to/B_rep1.fastq | 1     |
| B         | B_2       | /path/to/B_rep2.fastq | 1     |
| B         | B_3       | /path/to/B_rep3.fastq | 1     |

`order` column define which condition to treat as a control versus which treatment. `order = x` will be used as control for any `order > x`. In this case: condition B (treatment) versus A (control).

##### Library file

This file have to be comma-separeted.

First column is for guide's ID. This column should not be redundant.

Second column is guide's sequence.

Third column is the element targeted by the corresponding guide. Note: Multiple guides can target the same element.

```
guide_A.1,CTTAGTTTTGAACAAGTACA,element_A
guide_A.2,GTTGAGTTATCACACATCAT,element_A
guide_A.3,AATGTAGTGTAGCTACAGTG,element_A
guide_B.1,TTAGTTTATATCTTATGGCA,element_B
guide_B.2,GATTGTCTGTGAAATTTCTG,element_B
```

#### Starting from aligned BAM files

##### Design matrix

Let say that you have two conditions, A (control) and B (treatment), with 3 replicates for each but aligned BAM files instead of raw FASTQ files:

- Replace column name `fastq` by `bam` and fastq files paths by bam files paths, as follow:

| condition | replicate | bam                 | order |
|-----------|-----------|---------------------|-------|
| A         | A_1       | /path/to/A_rep1.bam | 0     |
| A         | A_2       | /path/to/A_rep2.bam | 0     |
| A         | A_3       | /path/to/A_rep3.bam | 0     |
| B         | B_1       | /path/to/B_rep1.bam | 1     |
| B         | B_2       | /path/to/B_rep2.bam | 1     |
| B         | B_3       | /path/to/B_rep3.bam | 1     |




#### Starting from count matrix

When starting from a count matrix, fastq/bam column isn't necessary:

| condition | replicate | order |
|-----------|-----------|-------|
| A         | A_1       | 0     |
| A         | A_2       | 0     |
| A         | A_3       | 0     |
| B         | B_1       | 1     |
| B         | B_2       | 1     |
| B         | B_3       | 1     |

However, its mandatory that `replicate` column contain the same labels that in count matrix header:

| shRNA       | Gene      | A_1 | A_2 | A_3 | B_1 | B_2 | B_3 |
|-------------|-----------|-----|-----|-----|-----|-----|-----|
| guide_A.1   | element_A | 456 | 273 | 345 | 354 | 587 | 258 |
| guide_A.2   | element_A | 354 | 234 | 852 | 546 | 64  | 452 |
