<img src="PitViper/docs/logo/pitviper_remasteredv2.png" alt="alt text" width="500" height="175">

> **An pipeline to analyze functional screening data from shRNA, CRISPR/Cas9 and CRISPR/dCas9 experiments.**

## Introduction

PitViper is intended to facilitate the analysis of functional screening data from various experiments (shRNA, CRISPR/Cas9 or CRISPR/dCas9).

The pipeline is built with [`Snakemake`](https://snakemake.readthedocs.io/en/stable/) and [`Flask`](https://flask.palletsprojects.com/en/2.0.x/), a workflow management system to create reproducible and scalable data analyses and a web framework.

### Prerequisites

To retrieve, install and run PitViper, we need [`Conda`](https://docs.conda.io/en/latest/) and [`Git`](https://git-scm.com/) available from the commande-line.

Once Conda is install and in order to speed-up installation process, we recommand to install [`Mamba`](https://github.com/mamba-org/mamba):

```bash
$  conda install -c conda-forge mamba
```

Now, Mamba should be available from the command-line.

Choose in wich directory you want to install PitViper and change your working directory for it and clone PitViper repository using `Git`:

```bash
$ mkdir myDir  # Create ./myDir/ directory in working directory

$ cd myDir  # Change working directory for ./myDir/

$ git clone https://github.com/PaulArthurM/PitViper.git  # Clone PitViper reposity in ./myDir/PitViper/

$ cd PitViper/PitViper  # Your are now in PitViper root directory
```

