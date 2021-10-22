<img src="PitViper/docs/logo/pitviper_remasteredv2.png" alt="alt text" width="500" height="175">

> **An pipeline to analyze functional screening data from shRNA, CRISPR/Cas9 and CRISPR/dCas9 experiments.**

## Introduction

PitViper is intended to facilitate the analysis of functional screening data from various experiments (shRNA, CRISPR/Cas9 or CRISPR/dCas9).

The pipeline is built with [`Snakemake`](https://snakemake.readthedocs.io/en/stable/), a workflow management system to create reproducible and scalable data analysis and [`Flask`](https://flask.palletsprojects.com/en/2.0.x/), a lightweight web framework.


##### Table of contents

- [Prerequisites](#prerequisites)

- [Set-up](#set-up)

- [Inputs](#inputs)


### Prerequisites

To retrieve, install and run PitViper, we need [`Conda`](https://docs.conda.io/en/latest/) and [`Git`](https://git-scm.com/) available from the commande-line.

Once Conda is installed and in order to speed-up installation process, we recommand to install [`Mamba`](https://github.com/mamba-org/mamba):

```bash
$  conda install -c conda-forge mamba
```

Now, Mamba should be available from the command-line.

Choose in wich directory you want to install PitViper and change your working directory for it, then clone PitViper repository at this location using `Git`:

```bash
$ mkdir myDir  # Create ./myDir/ directory in working directory

$ cd myDir  # Change working directory for ./myDir/

$ git clone https://github.com/PaulArthurM/PitViper.git  # Clone PitViper reposity in ./myDir/PitViper/

$ cd PitViper/PitViper  # Your are now in PitViper root directory
```

### Set-up

Then, we need to install PitViper dependancies. To facilitate dependancy management, a Conda YAML file containing all dependancies have been created and can be used along with `Mamba` to automatically install them: `pitviper_env.yaml`.

Furthermore, `install_PitViper_env.sh` is a bash script created to perform this step in one command:

```bash
$ ./install_PitViper_env.sh
```

`./install_PitViper_env.sh` first create a Conda environment call `pitviper_env` that should be visible as follow if you run:

```bash
$ conda env list
  # conda environments:
  #
  base                  *  /home/paularthur/miniconda3
  pitviper_env             /home/paularthur/miniconda3/envs/pitviper_env
```

Once `pitviper_env` is created, you can run `run.sh` script, it should open PitViper GUI in your default web browser.

### Inputs

Prior to use PitViper, we need to introduce it's input files.

PitViper allow to start an analysis from raw data files such as FASTQ, already aligned BAM file or a count matrix. Depending of that, input files will be different.

#### Starting from raw FASTQ files

You will need:

1. Path to FATSQ files on system
2. A library file with three comma-separated columns without header: shRNA ID, shRNA sequence, targer element.
3. A design matrix that summary your conditions, their replicates and associated FASTQ files.

##### Design matrix

Let say that you have two conditions, A and B, with 3 replicates for each. The associated matrix should look as below:

| condition | replicate | fastq                 |
|-----------|-----------|-----------------------|
| A         | A_1       | /path/to/A_rep1.fastq |
| A         | A_2       | /path/to/A_rep2.fastq |
| A         | A_3       | /path/to/A_rep3.fastq |
| B         | B_1       | /path/to/B_rep1.fastq |
| B         | B_2       | /path/to/B_rep2.fastq |
| B         | B_3       | /path/to/B_rep3.fastq |

Those name will be used for differentiation in results report. 

##### Library file

First column indicate the ID of the shRNA. This ID in this column should not be redundant.

Second column is shRNA sequence.

Thirs column indicate the element targeted by the corresponding shRNA sequence. Note: Multiple shRNA can target the same element.

```
element_A.1,CTTAGTTTTGAACAAGTACA,element_A
element_A.2,GTTGAGTTATCACACATCAT,element_A
element_A.3,AATGTAGTGTAGCTACAGTG,element_A
element_B.1,TTAGTTTATATCTTATGGCA,element_B
element_B.2,GATTGTCTGTGAAATTTCTG,element_B
```

