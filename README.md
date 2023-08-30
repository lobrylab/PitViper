<img src="PitViper/docs/logo/pitviper_remasteredv2.png" alt="alt text" width="500" height="175">

> **A pipeline to analyze functional screening data from shRNA, CRISPR/Cas9 and CRISPR/dCas9 experiments.**

## Introduction
Processing, InTerpretation and VIsualization of PoolEd screening Results

PitViper is intended to facilitate the analysis of functional screening data from various functionnal screening experiments (shRNA, CRISPR/Cas9 or CRISPR/dCas9).

The pipeline is built with [`Snakemake`](https://snakemake.readthedocs.io/en/stable/), a workflow management system to create reproducible and scalable data analysis, [`Flask`](https://flask.palletsprojects.com/en/2.0.x/), a lightweight web framework and [`Jupyter`](https://jupyter.org/), a web application for creating and sharing computational documents.


##### Table of contents

- [Prerequisites](#prerequisites)

- [Installation](#installation)

  - [From the command line](#from-the-command-line)
  
  - [Run PitViper from Docker container](#run-pitviper-from-docker-container)

- [Inputs](#inputs)

  - [Starting from raw FASTQ or BAM files](#starting-from-raw-fastq-files)

  - [Starting from count matrix](#starting-from-count-matrix)


### Prerequisites

[`Conda`](https://docs.conda.io/en/latest/) and [`Git`](https://git-scm.com/) are required to install PitViper.

Then, clone PitViper repository using `Git`:

```bash
$ git clone https://github.com/PaulArthurM/PitViper.git  # Clone PitViper reposity in ~/PitViper/

$ cd PitViper  # Change working directory to PitViper root directory
```

### Installation

#### From the command line

The `run.sh` script can be used to automatically create the conda environment, install all dependencies and run PitViper graphical user interface (GUI). Dependencies are installed only when `run.sh` is run for the first time.

```bash
# Inside PitViper root directory
$ ./run.sh
```

After installation is completed, the PitViper GUI should be automatically opened in your default web browser. 

<img src="PitViper/docs/PitViper.png" alt="alt text">

#### Run PitViper from Docker container

##### Building Image

PitViper main directory contains a Dockerfile that allows easy building of PitViper image. From this folder simply run:

```bash
$ docker build -t [image_name] .
```

##### Downloading image from dockerhub

PitViper docker image can also be downloaded from dockerhub using:

```bash
$ docker pull lobrylab/pitviper:v1.0
```

##### Running container

To start a PitViper container simply run:

```bash
$ docker run -p 5000:5000 -p 8888:8888 -v [fastq/bam_files_path]: [fastq/bam_files_path] -n [name_the_container] [PitViper_image_name]
```

For example:

```bash
$ docker run -p 5000:5000 -p 8888:8888 -v /home/data/:/home/data/ -n pitviper lobrylab/pitviper:v1.0
```

PitViper GUI can now be accessed in your web browser at the address: localhost:5000.

Upon completion of PitViper analysis, the jupyter notebook will be accessible following the localhost:8888/[token] link displayed in the terminal.

To quit PitViper and stop the container simply quit the jupyter notebook session.

##### Re-starting container

To start a new analysis, just restart the container using the following command :

```bash
$ docker start -a [container_name]
```

##### Accessing jupyter notebooks of a previous docker PitViper runs

To access the notebook of a previously ran PitViper analysis, first start the container :

```bash
$ docker start [container_name]
```

Then execute the notebook script: 

```bash
$ docker exec -ti [container_name] bash /PitViper/notebook.sh
```

The folder containing all previously ran analysis will be accessible in jupyter notebook following the localhost:8888/[token] link displayed in the terminal.

### Inputs

Before using PitViper, we will introduce its input files.

PitViper allows you to start an analysis from raw data files such as FASTQ, an already aligned BAM file or a count matrix. Depending on this, the input files will be different.

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
