# SynDESeq2

This repository contains snakemake workflows used to perform differential gene expression analysis using DESeq2 and 
GO Term analysis via ClusterProfiler. 

## Hardware requirements
SynDESeq2 only requires a standard computer with sufficient RAM depending on the analyzed dataset.

## Software requirements

### OS
SynDESeq2 is supported for Linux. 

### Dependencies

The software was tested on Ubuntu 22.04 using following dependencies:

```text
conda
python=3.12.10,
snakemake
```


## Install

Fetch the Git repository and install the following dependencies

```shell
conda create -n SynDESeq2
conda activate SynDESeq2
conda install python==3.12 snakemake
```

## Running

### Demo
You can run the pipeline using the axample `airway` dataset via the following command:

```shell
snakemake -s snakefile.smk  --cores 1 --configfile example_config.yml --use-conda --conda-frontend conda
```

Note that you need to download the Git repository for that


### Using your data

You can either create a custom config.yaml file and run the Workflow as shown in the Demo or include the workflow from 
github in your own snakemake file as follows:

```python
snakefile = github("domonik/SynDESeq2", path="snakefile.smk")

module syndeseq_workflow:
    snakefile: snakefile
    config: config["SynDESeqConfig"]

use rule * from syndeseq_workflow as syndeseq_*
```

This way you have access to all the rules in this workflow, or you can just generate the files via requesting them in 
 your target rule:


```python
rule all:
    default_target: True
    input:
        syndeseq = rules.syndeseq_all.input,
```


## Visualization

The results of a workflow run can be visualized using the DEplots python package. Visit 
 https://github.com/domonik/DEplots for more details.


