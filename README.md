# 2025-pi-disco-pub

[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)


## Purpose

This repo accompanies the pub ["Case study: Using AlphaFold-Multimer to predict the targets of tick protease inhibitors"](https://doi.org/10.57844/arcadia-77d4-1c5d). We are interested in developing approaches to predict the targets of tick secreted effector proteins, with the goal of learning new strategies for modulating itch, pain, and inflammation in the skin. In this work, we try out using AlphaFold Multimer to predict targets for 10 representatives of a family of tick TIL-domain protease inhibitors (orthogroup OG0000058). 

This pilot had three main steps:
1. **Protease inhibitor discovery**: Identify protease inhibitor gene families across 15 species of ticks using [ProteinCartography](https://github.com/Arcadia-Science/ProteinCartography) and orthogroup data from [previous work](https://dx.doi.org/10.57844/arcadia-4e3b-bbea)
2. **Phylogenomic profiling**: Identify which gene family is most predictive of the ability of ticks to suppress host-detection mechanisms such as itch, pain, and inflammation using phylogenomic profiling 
3. **Protein-protein interaction prediction**: Predict targets for protease inhibitor family orthogroup OG0000058 using AlphaFold Multimer

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following commands to create the pipeline run environments.

To create the enviroment used for protease inhibitor discovery:

```{bash}
mamba env create -n pi-disco --file envs/pi-disco.yml
conda activate pi-disco
```
To create the enviroment used for protein-protein interaction prediction:

```{bash}
mamba env create -n tick_ppi --file envs/ppi.yml
conda activate tick_ppi
```

Our phylogenomic profiling work ([Notebook 04](https://github.com/Arcadia-Science/2025-pi-disco-pub/blob/main/notebooks/04_trait_mapping.Rmd)) uses R scripts and an R notebook. The R scripts can either be run directly in the terminal as R scripts or in the terminal R environment, or can be run in RStudio. If you choose to run these scripts in Rstudio, you can install it by following the instructions [here](https://posit.co/download/rstudio-desktop/). Once RStudio is installed, run the following from command line to open Rstudio from your activated conda environment:

```{bash}
open -a Rstudio
```

To execute the code in [Notebook 04](https://github.com/Arcadia-Science/2025-pi-disco-pub/blob/main/notebooks/04_trait_mapping.Rmd), please make sure you have installed all the dependencies shown in the ```sessionInfo()``` at the beginning the [session](notebooks/04_trait_mapping.nb.html)

## Data related to protease inhibitor discovery  
- Datasheets used in Notebooks 01-03 and 05 are found in [datasheets](https://github.com/Arcadia-Science/2025-pi-disco-pub/tree/main/datasheets)
- All chelicerate protein annotations are in [all_chelicerate_proteins.csv](https://zenodo.org/records/15186244) on Zenodo
- Tick protein structures are in [pi-disco_esmfold_structures.zip](https://zenodo.org/records/15186244)
- Animal toxin database downloaded from [UniProt Animal toxin annotation project](https://www.uniprot.org/help/Toxins) 
- Our full ProteinCartography run is in [tick_PUFs_PIs_1200_plus_toxinDB_carto_run3.zip](https://zenodo.org/records/15186244) on Zenodo, key results are found in [outputs/carto_run](https://github.com/Arcadia-Science/2025-pi-disco-pub/tree/main/outputs/carto_run)
- [NovelTree orthogroup data](https://github.com/Arcadia-Science/2025-pi-disco-pub/blob/main/datasheets/Orthogroups.tsv) comes from previous work [here](https://dx.doi.org/10.57844/arcadia-4e3b-bbea)

## Data related to phylogenomic profiling
- Datasheets of [high-confidence](https://github.com/Arcadia-Science/2025-pi-disco-pub/blob/main/datasheets/PI_orthogroups_high_qual_06112023.tsv) and [low-confidence protease](https://github.com/Arcadia-Science/2025-pi-disco-pub/blob/main/datasheets/PI_orthogroups_low_qual_06112023.tsv) inhibitor gene families
- NovelTree workflow output directory [chelicerata-v1-10062023](https://zenodo.org/records/14113178)
- Results from phylogenomic profiling in [detection_suppression_outputs](https://github.com/Arcadia-Science/2025-pi-disco-pub/tree/main/outputs/detection_suppression_outputs/trait_prediction/detection_suppression_test_lasso/s_counts_proteases_combined)

## Data related to Protein-protein interaction prediction
- Datasheets of Uniprot accessions and expression metadata found in [datasheets](https://github.com/Arcadia-Science/2025-pi-disco-pub/tree/main/datasheets)
- Results from D-script and AF-multimer runs are in [protein_protein_interaction_results](https://github.com/Arcadia-Science/2025-pi-disco-pub/tree/main/outputs/protein_protein_interaction_results)
- Raw results files from AF-Multimer are on [Zenodo](https://zenodo.org/records/15186244) The .json and .pae files are the most important for reanalyzing the data to calculate AF-Multimer metrics

## Overview

### Description of the folder structure
- `datasheets/` contains input metadata and datasheets used in analysis notebooks 
- ` envs/` contains two ```.yml``` files that specify the conda environments used for protease inhibitor discovery and protein-protein interaction prediction
- `notebooks/` contains the Jupyter notebooks used to do this analysis 
- `outputs/` contains figures, select results from the ProteinCartography run, results from the phylogenomic profiling pipeline, and results from protein-protein interaction prediction
- `scripts/` contains: scripts for folding proteins using an Arcadia-hosted ESMFold API, further details in `scripts/README_esmfold_scripts.md`, scripts for protein-protein interaction prediction using D-SCRIPT and AlphaFold multimer, with further details in `scripts/README_ppi_scripts.md`, and an R script `scripts/trait_mapping_functions.R` that accompanies [Notebook 04](https://github.com/Arcadia-Science/2025-pi-disco-pub/blob/main/notebooks/04_trait_mapping.Rmd). 
  

### Methods

Protease inhibitor discovery:
1. TO DO

Phylogenomic profiling:
1. TO DO

Protein-protein interaction prediction: 
1. Activate environment `envs/ppi.yml`
2. Collect sets of tick and human proteins. Human proteins and associated metadata are pulled down from Uniprot accessions with the [`fetch_accession.py`](https://github.com/Arcadia-Science/ProteinCartography/blob/v0.4.0-alpha/ProteinCartography/fetch_accession.py) and [`fetch_uniprot_metadata.py`](https://github.com/Arcadia-Science/ProteinCartography/blob/v0.4.0-alpha/ProteinCartography/fetch_uniprot_metadata.py) from [ProteinCartography](https://github.com/Arcadia-Science/ProteinCartography) following the instructions there and providing TXT files of Uniprot accessions.
3. Make sequence-based predictions of PPIs with D-script using the script `make_dscript_PPI_predictions.py`. Analysis of these results and attempted filtering are in `scripts/dscript-results-analysis.R`.
4. Prepare FASTA files for running with the Goolge Colab Batch AF-multimer notebook. This can be done either providing the separate bait and prey FASTA sequences that should be combined, or from a TSV of comparisons and a FASTA containing all the sequences
5. Analyze Alphafold results with the Jupyter notebook code in [`notebooks/06_calculate_afmultimer_predictions.ipynb`](./notebooks/06_calculate_afmultimer_predictions.ipynb) and run with the [LazyAF Google Colab notebook](https://colab.research.google.com/drive/1j7WJLcUHTR8BrjkWDaU549rFk6X5Zu18), but paste in the modified code in this notebook. The modified code adds lines to calculate the average pLDDT and include in the resulting CSV files for each comparison. Change the text in each cell to point to your data and rename the output CSV file. If you don't need the pLDDT information just use the default notebook code.
6. Save both the results and analysis folders from the AF-multimer predictions and the LazyAF analysis, which selects the highest ranking prediction JSON file and copies it into the analysis folder. Downloading the folder from Google Drive will automatically create a .zip archive and save this in your results folder for backup.
7. Analyze AF-multimer results compared to D-script predictions with [`notebooks/07_analyze-dscript-afmultimer-comparisons.ipynb`](./notebooks/07_analyze-dscript-afmultimer-comparisons.ipynb)
8. Code for analysis and plotting is in the [`notebooks/08_full-filtered-serine-protease-analysis.ipynb`](./notebooks/08_full-filtered-serine-protease-analysis.ipynb)
9. We calculated additional AFmultimer metrics with the [`colabfold_analysis.py`](https://github.com/walterlab-HMS/AF2multimer-analysis/blob/main/colabfold_analysis.py) script as part of the AF2multimer-analysis toolkit that underlies the [predictomes.org](https://predictomes.org/) webserver. To use the script:
   1.  Download the raw AF2multimer results from Google drive and unzip them. Google will sometimes split up large downloads into multiple zip files and not all results files that go together are in the same subdirectories. Make sure for this script to work all the `.json` and `.pae` files together as well as a `.done.txt` file with the name of the comparison
   2.  Clone the directory with `git clone https://github.com/walterlab-HMS/AF2multimer-analysis`
   3.  Run the script with `python3 colabfold_analysis.py <AF2MULTIMER_RESULTS_DIRECTORY>`. This will create a new folder with the same name as your input directory with `_analysis` appended at the end. This creates three CSV files: a contacts, interfaces, and summary file. The `summary.csv` is the most useful for downstream steps and is how we calculated the pDockQ score for each interaction prediction. The contacts and interfaces files give more information for each residue the calculated metrics, whereas the summary file is an average across all residues.

### Compute Specifications

Scripts and analysis run on MacOS except for AF-multimer runs which were through Google Colab. Most of the runs on Google Colab were with A100 GPUs, and when not available with either V100 or T4 GPUs.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md). 
