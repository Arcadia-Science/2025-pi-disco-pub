# 2025-pi-disco-pub

[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)


## Purpose

TODO: Briefly describe the core analyses performed in the repository and the motivation behind them.

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
TODO: Replace <NAME> with the name of your environment
mamba env create -n <NAME> --file envs/dev.yml
conda activate <NAME>
```

## Data

TODO: Add details about the description of input / output data and links to Zenodo depositions, if applicable.

## Overview

### Description of the folder structure

### Methods

TODO: Include a brief, step-wise overview of analyses performed.

> Example:
>
> 1.  Download scripts using `download.ipynb`.
> 2.  Preprocess using `./preprocessing.sh -a data/`
> 3.  Run Snakemake pipeline `snakemake --snakefile Snakefile`
> 4.  Generate figures using `pub/make_figures.ipynb`.

### Compute Specifications

TODO: Describe what compute resources were used to run the analysis. For example, you could list the operating system, number of cores, RAM, and storage space.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).




# Protease inhibitor protein-protein interactions

## Purpose

This repository is for analyzing predictions of tick-human protein-protein interactions (PPIs). Documentation for how to use the command-line scripts is provided in the README of the `scripts` directory. Other scripts are for analyzing or organizing metadata and don't have command-line documentation but brief descriptions of what they were used for. Analysis is primarily in different notebooks in the `notebooks` subdirectory.

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

To use the scripts in this repository create a conda environment:
```
mamba env create -n tick_ppi --file envs/ppi.yml
mamba activate tick_ppi
```

Additionally this project used a pip package D-script to predict sequence-based PPI interactions. This tool is a structure informed, deep-learning method trained on experimentally verified protein-protein interactions and gives confidence scores to combinations of proteins for predicting if they might interact based on sequence alone. To make PPI predictions with this tool you will need to download a model. For this project we used the recommended Human Topsy-Turvy model available [here](http://cb.csail.mit.edu/cb/dscript/data/models/topsy_turvy_v1.sav). This model is of human-human protein-protein interactions and was trained as a sequence-based, deep-learning model for PPI prediction. Read more about the Topsy-Turvy model [here](https://cb.csail.mit.edu/cb/topsyturvy/). There are other models for other within-organism protein-protein interactomes, such as yeast, mice, etc.

## Running AF-multimer in Google Colab

Using either the [`make-afmultimer-fasta-from-tsv.py`](./scripts/make-afmultimer-fasta-from-tsv.py) or [`make-afmultimer-fasta-from-seqs.py`](./scripts/make-afmultimer-fasta-from-seqs.py) to make the FASTA files compatible for input to Alphafold multimer, upload these FASTA files to a folder in Google Drive. Run AF-multimer Batch in Google Colab [here](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/batch/AlphaFold2_batch.ipynb). The only modifications to the existing code are mounting Google Drive, changing the paths in cell 1 from `drive.mount('/content/drive')` to `drive.mount('/content/gdrive')`. Accordingly change this in cell 2 to point the path for where your folder of FASTA sequences are, such as `/content/gdrive/MyDrive/serine_proteases_under_700`. When you mount to your Google Drive (this can be any google drive folder) you will grant read/write permissions either to the entire contents of your entire Google Drive or select folders that you are reading/writing to with the AF-multimer Colab.

## Data

Input data are sets of tick and/or human serine proteases from Uniprot. Human proteins from Uniprot were pulled down for accessions listed in the `metadata` Uniprot table. Tick proteins were either selected from the [host detection suppression trait association analysis](https://github.com/Arcadia-Science/2024-chelicerate-phylogenomics) or from the analysis described above. These combinations of proteins were used for both sequence-based and structure-based predictions of PPIs. 

### Description of the folder structure
- `metadata/` Lists of Uniprot accessions and descriptions of annotations
- `notebooks/` Analysis notebook code used in Google Colab to collate results from AF-multimer, and notebooks for comparing results from D-script and AF-multimer and analyzing and plotting protease inhibitor results
- `scripts/` Command-line python scripts for generating AF-multimer FASTAs and automating D-script predictions, and analysis scripts
- `figs/` Figures produced
- `results/` Results, including files from D-script and AF-multimer

The `results/` directory has results from D-SCRIPT and AF-Multimer runs. The final filtered results table with metrics from AF-Multimer, annotation data, and expression data is in [`serine_protease_results/2024-04-24-final-serine-protease-filtered-interface-comps-results.tsv`](./results/serine_protease_results/2024-04-24-final-serine-protease-filtered-interface-comps-results.tsv). The files here are used in the analysis notebooks and also described in the notebooks.

Raw results files from AF-Multimer are stored in an S3 Bucket on AWS since these are large files: https://us-west-1.console.aws.amazon.com/s3/buckets/2024-tick-ppi?region=us-west-1. The `.json` and `.pae` files are the most important for reanalyzing the data to calculate AF-Multimer metrics.

### Methods

1. Collect sets of tick and human proteins. Human proteins and associated metadata are pulled down from Uniprot accessions with the [`fetch_accession.py`](https://github.com/Arcadia-Science/ProteinCartography/blob/v0.4.0-alpha/ProteinCartography/fetch_accession.py) and [`fetch_uniprot_metadata.py`](https://github.com/Arcadia-Science/ProteinCartography/blob/v0.4.0-alpha/ProteinCartography/fetch_uniprot_metadata.py) from [ProteinCartography](https://github.com/Arcadia-Science/ProteinCartography) following the instructions there and providing TXT files of Uniprot accessions.
2. Make sequence-based predictions of PPIs with D-script using the script `make_dscript_PPI_predictions.py`. Analysis of these results and attempted filtering are in `scripts/dscript-results-analysis.R`.
3. Prepare FASTA files for running with the Goolge Colab Batch AF-multimer notebook. This can be done either providing the separate bait and prey FASTA sequences that should be combined, or from a TSV of comparisons and a FASTA containing all the sequences
4. Analyze Alphafold results with the Jupyter notebook code in [`notebooks/calculate_afmultimer_predictions.ipynb`](./notebooks/calculate_afmultimer_predictions.ipynb) and run with the [LazyAF Google Colab notebook](https://colab.research.google.com/drive/1j7WJLcUHTR8BrjkWDaU549rFk6X5Zu18), but paste in the modified code in this notebook. The modified code adds lines to calculate the average pLDDT and include in the resulting CSV files for each comparison. Change the text in each cell to point to your data and rename the output CSV file. If you don't need the pLDDT information just use the default notebook code.
5. Save both the results and analysis folders from the AF-multimer predictions and the LazyAF analysis, which selects the highest ranking prediction JSON file and copies it into the analysis folder. Downloading the folder from Google Drive will automatically create a .zip archive and save this in your results folder for backup.
6. Analyze AF-multimer results compared to D-script predictions with [`notebooks/analyze-dscript-afmultimer-comparisons.ipynb`](./notebooks/analyze-dscript-afmultimer-comparisons.ipynb)
7. Code for analysis and plotting is in the [`notebooks/full-filtered-serine-protease-analysis.ipynb`](./notebooks/full-filtered-serine-protease-analysis.ipynb)
8. Calculate all-v-all TMscore comparisons with `foldseek`, where we did this for comparisons of tick protease inhibitor and human serine protease structural similarity, respectively. This was done with the `foldseek easy-search` subworkflow
9. We calculated additional AFmultimer metrics with the [`colabfold_analysis.py`](https://github.com/walterlab-HMS/AF2multimer-analysis/blob/main/colabfold_analysis.py) script as part of the AF2multimer-analysis toolkit that underlies the [predictomes.org](https://predictomes.org/) webserver. To use the script:
   1.  Download the raw AF2multimer results from Google drive and unzip them. Google will sometimes split up large downloads into multiple zip files and not all results files that go together are in the same subdirectories. Make sure for this script to work all the `.json` and `.pae` files together as well as a `.done.txt` file with the name of the comparison
   2.  Clone the directory with `git clone https://github.com/walterlab-HMS/AF2multimer-analysis`
   3.  Run the script with `python3 colabfold_analysis.py <AF2MULTIMER_RESULTS_DIRECTORY>`. This will create a new folder with the same name as your input directory with `_analysis` appended at the end. This creates three CSV files: a contacts, interfaces, and summary file. The `summary.csv` is the most useful for downstream steps and is how we calculated the pDockQ score for each interaction prediction. The contacts and interfaces files give more information for each residue the calculated metrics, whereas the summary file is an average across all residues.

### Note

The majority of serine protease predictions were made using a single **_Amblyomma_ serine protease inhibitor**. The serine proteases/tick inhibitor runs produced about 50 high-confidence predictions (includes reviewed and unreviewed Uniprot hits), and there is [a preprint where analysis of cross-kingdom plant-bacterial pathogen protease inhibitor interactions produces several high-confidence predictions](https://www.biorxiv.org/content/10.1101/2023.04.03.535425v1.full). The high volume of high-confidence predictions may be because there are experimental crystal structures of protease inhibitor complexes, specifically for plant-bacterial pathogen pairs.

### Compute Specifications

Scripts and analysis run on MacOS except for AF-multimer runs which were through Google Colab. Most of the runs on Google Colab were with A100 GPUs, and when not available with either V100 or T4 GPUs.
