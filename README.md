# 2025-pi-disco-pub

[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)


## Purpose

This repo accompanies the pub ["Case study: Using AlphaFold-Multimer to predict the targets of tick protease inhibitors"](https://doi.org/10.57844/arcadia-77d4-1c5d). We are interested in developing approaches to predict the targets of tick secreted effector proteins, with the goal of learning new strategies for modulating itch, pain, and inflammation in the skin. In this work, we try out using AlphaFold Multimer to predict targets for 10 representatives of a family of tick TIL-domain protease inhibitors (orthogroup OG0000058). 

This pilot had three main steps:
1. **Protease inhibitor discovery**: Identify protease inhibitor gene families across 15 species of ticks using [ProteinCartography](https://github.com/Arcadia-Science/ProteinCartography) and orthogroup data from [previous work](https://dx.doi.org/10.57844/arcadia-4e3b-bbea)
2. **Phylogenomic profiling**: Identify which gene family is most predictive of the ability of ticks to suppress host-detection mechanisms such as itch, pain, and inflammation using trait mapping
3. **Protein-protein interaction prediction**: Predict targets for protease inhibitor family orthogroup OG0000058 using AlphaFold Multimer. We also tried D-SCRIPT, but ultimately didn't move forward with that approach.

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

Our phylogenomic profiling workflow is implemented in an R Markdown notebook ([Notebook 04](./notebooks/04_trait_mapping.Rmd)), which sources a ```trait_mapping_functions.R``` script located in the [`scripts`](./scripts) folder. To run this analysis, you can install Rstudio by following the instructions [here](https://posit.co/download/rstudio-desktop/). Once RStudio is installed, run the following from command line to open Rstudio from your activated conda environment:

```{bash}
open -a Rstudio
```

To execute the code in [Notebook 04](./notebooks/04_trait_mapping.Rmd), please make sure you have installed all the dependencies shown in the ```sessionInfo()``` printed at the beginning the [notebook](./notebooks/04_trait_mapping.nb.html).

## Data related to protease inhibitor discovery  
- Datasheets used in Notebooks 01-03 and 05 are found in [`datasheets`](./datasheets)
- All chelicerate protein annotations are in [`all_chelicerate_proteins.csv`](https://zenodo.org/records/15186244) on Zenodo
- Tick protein structures are in [`pi-disco_esmfold_structures.zip`](https://zenodo.org/records/15186244)
- Animal toxin database downloaded from [UniProt Animal toxin annotation project](https://www.uniprot.org/help/Toxins) 
- Our full ProteinCartography run is in [`tick_PUFs_PIs_1200_plus_toxinDB_carto_run3.zip`](https://zenodo.org/records/15186244) on Zenodo, key results are found in [`outputs/carto_run`](./outputs/carto_run)
- [NovelTree orthogroup data](./datasheets/Orthogroups.tsv) comes from previous work [here](https://dx.doi.org/10.57844/arcadia-4e3b-bbea)

## Data related to phylogenomic profiling
- Datasheets of [high-confidence](./datasheets/PI_orthogroups_high_qual_06112023.tsv) and [low-confidence protease](./datasheets/PI_orthogroups_low_qual_06112023.tsv) inhibitor gene families
- NovelTree workflow output directory [`chelicerata-v1-10062023`](https://zenodo.org/records/14113178)
- Results from phylogenomic profiling in [`detection_suppression_outputs`](./outputs/detection_suppression_outputs/trait_prediction/detection_suppression_test_lasso/s_counts_proteases_combined)

## Data related to Protein-protein interaction prediction
- Datasheets of Uniprot accessions and expression metadata found in [`datasheets`](./datasheets)
- Results from D-script and AF-multimer runs are in [`protein_protein_interaction_results`](./outputs/protein_protein_interaction_results)
- Raw results files from AF-Multimer are on [Zenodo](https://zenodo.org/records/15186244) The .json and .pae files are the most important for reanalyzing the data to calculate AF-Multimer metrics

## Overview

### Description of the folder structure
- [`datasheets/`](./datasheets) contains input metadata and datasheets used in analysis notebooks 
- [`envs/`](./envs) contains two `.yml` files that specify the conda environments used for protease inhibitor discovery and protein-protein interaction prediction
- [`notebooks/`](./notebooks) contains the Jupyter notebooks used to do this analysis 
- [`outputs/`](./outputs) contains [figures](./outputs/figures), [select results](./outputs/carto_run) from the ProteinCartography run, [results](./outputs/detection_suppression_outputs/trait_prediction/detection_suppression_test_lasso/s_counts_proteases_combined) from the phylogenomic profiling pipeline, and [results](./outputs/protein_protein_interaction_results) from protein-protein interaction prediction
- [`scripts/`](./scripts) contains: scripts for folding proteins using an Arcadia-hosted ESMFold API, further details in [`scripts/README_esmfold_scripts.md`](./scripts/README_esmfold_scripts.md), scripts for protein-protein interaction prediction using D-SCRIPT and AlphaFold multimer, with further details in [`scripts/README_ppi_scripts.md`](./scripts/README_ppi_scripts.md), and an R script [`scripts/trait_mapping_functions.R`](./scripts/trait_mapping_functions.R) that accompanies [Notebook 04](./notebooks/04_trait_mapping.Rmd). 
  

### Methods

**Protease inhibitor discovery:**
1. Activate environment [`envs/pi-disco.yml`](./envs/pi-disco.yml).
2. Get proteins: [`01_get_proteins.ipynb`](./notebooks/01_get_proteins.ipynb) \
Select protease inhibitors (PIs) and secreted proteins of unknown function (PUFs) <1200 aa, using annotations found in [`all_chelicerate_proteins.csv`](https://zenodo.org/records/15186244). PIs that we have selected are in [`./datasheets/tick_PIs_1200.csv`](./datasheets/tick_PIs_1200.csv) and secreted PUFs are in [`./datasheets/tick_PUFs_1200.csv`](./datasheets/tick_PUFs_1200.csv).
3. Fold tick proteins for structural analysis: 
We use an Arcadia-hosted ESMFold API to predict and download the folded proteins. All the scripts we used are in [`./scripts`](./scripts), and detailed usage can be found in [`scripts/README_esmfold_scripts.md`](./scripts/README_esmfold_scripts.md). The structures that we computed can be found in [`pi-disco_esmfold_structures.zip`](https://zenodo.org/records/15186244).
4. Download structures of Animal Toxins:
We downloaded structures from [UniProt Animal toxin annotation project](https://www.uniprot.org/help/Toxins) using [ProteinCartography](https://github.com/Arcadia-Science/ProteinCartography/tree/v0.4.0-alpha) scripts [`fetch_accession.py`](https://github.com/Arcadia-Science/ProteinCartography/blob/v0.4.0-alpha/ProteinCartography/fetch_accession.py) and [`fetch_uniprot_metadata.py`](https://github.com/Arcadia-Science/ProteinCartography/blob/v0.4.0-alpha/ProteinCartography/fetch_uniprot_metadata.py). Our list of successfully downloaded structures can be found [here](./datasheets/toxinDB_pdb_list.csv). 
5. Prepare metadata files to run ProteinCartography on tick proteins and other animal toxins: [`02_make_metadata.ipynb`](./notebooks/02_make_metadata.ipynb) \
This notebook prepares the [metadata file](./datasheets/tick_PUFs_PIs_1200_plus_toxinDB_metadata.tsv) for tick proteins and toxin database proteins to go into ProteinCartography. 
6. Use ProteinCartography to cluster tick proteins and toxins:
We cloned the [ProteinCartography repo](https://github.com/Arcadia-Science/ProteinCartography/tree/v0.4.0-alpha), and ran ProteinCartography on the combined tick and toxin database structures using `cluster` mode. Our config file can be found [here](./outputs/carto_run/config_ff_run3.yml) and our full results can be found in [`tick_PUFs_PIs_1200_plus_toxinDB_carto_run3.zip`](https://zenodo.org/records/15186244). While the full dataset is too large to host on GitHub, we've moved the key files into an output folder [here](./outputs/carto_run/carto_output_run3/clusteringresults). 
7. Analyze ProteinCartography results to get protease inhibitor gene families:[03_analyze_cartography.ipynb](./notebooks/03_analyze_cartography.ipynb) \
This notebook takes in the clustered output data from the ProteinCartography run on the tick proteins combined with the toxin database. By looking at semantic data (found in [`tick_PUFs_PIs_1200_plus_toxinDB_carto_run3.zip`](https://zenodo.org/records/15186244)) and intra-cluster structural similarity, we find 7 clusters that are high quality protein inhibitors (2.5k proteins) and 4 clusters that are interesting but are much lower confidence (1.8 k proteins). We then combine these protease inhibitor sturctural clusters with NovelTree [orthogroup data](./datasheets/Orthogroups.tsv) to select protease inhibitor orthogroups to put through trait mapping. The notebook selects 36 orthogroups that have either [high quality](./datasheets/PI_orthogroups_high_qual_06112023.tsv) or [low quality](./datasheets/PI_orthogroups_low_qual_06112023.tsv) protease inhibitor annotations for trait mapping. 

**Phylogenomic profiling:**
1. Download NovelTree workflow output directory [`chelicerata-v1-10062023`](https://zenodo.org/records/14113178).
2. In R-Studio or R Jupyter Notebooks, run [Notebook 04](./notebooks/04_trait_mapping.Rmd). This performs trait mapping on the 36 orthogroups that have either [high quality](./datasheets/PI_orthogroups_high_qual_06112023.tsv) or [low quality](./datasheets/PI_orthogroups_low_qual_06112023.tsv) protease inhibitor annotations, using evolutionary data from the NovelTree analysis. You can see our full session [here](./notebooks/04_trait_mapping.nb.html), and full results can be found [here](./outputs/detection_suppression_outputs/trait_prediction/detection_suppression_test_lasso/s_counts_proteases_combined).
3. To dig further into the results and select proteins for PPI prediction, activate the [`pi-disco` environment](./envs/pi-disco.yml). Launch [Notebook 05](./notebooks/05_candidate_selection.ipynb), where we use model variable importance and individual conditional expectation (ICE) to select the gene families that most strongly and postively predict the detection-supression trait. We also filter by secretion signal and by expression level in the female _Amblyomma americanum_ salivary gland to select the most promsing individual proteins to study using protein-protein interaction predicition. We select 10 TIL domain proteins, which can be found [here](./datasheets/Amblyomma-americanum_OG0000058_candidates.tsv).

**Protein-protein interaction prediction:** 
1. Activate environment [`envs/ppi.yml`](./envs/ppi.yml)
2. Collect sets of tick and human proteins. Human proteins and associated metadata are pulled down from Uniprot accessions with the [`fetch_accession.py`](https://github.com/Arcadia-Science/ProteinCartography/blob/v0.4.0-alpha/ProteinCartography/fetch_accession.py) and [`fetch_uniprot_metadata.py`](https://github.com/Arcadia-Science/ProteinCartography/blob/v0.4.0-alpha/ProteinCartography/fetch_uniprot_metadata.py) from [ProteinCartography](https://github.com/Arcadia-Science/ProteinCartography) following the instructions there and providing TXT files of Uniprot accessions. All tick protein sequences used in this study can be found in [`2024-06-24-all-chelicerate-noveltree-proteins.fasta`](https://zenodo.org/records/14113178). We specifically looked at these 10 proteins [here](./datasheets/Amblyomma-americanum_OG0000058_candidates.tsv), using Amblyomma-americanum_evm.model.contig-94090-1.4 as a representative, since it had one of the highest-quality structures of the group (pLDDT = 83.8).
3. Make sequence-based predictions of PPIs with D-script using the script [`make_dscript_PPI_predictions.py`](./scripts/make_dscript_PPI_predictions.py). Analysis of these results and attempted filtering are in [`scripts/dscript-results-analysis.R`](scripts/dscript-results-analysis.R).
4. Prepare FASTA files for running with the Goolge Colab Batch AF-multimer notebook. This can be done either providing the separate bait and prey FASTA sequences that should be combined, or from a TSV of comparisons and a FASTA containing all the sequences
5. Analyze Alphafold results with the Jupyter notebook code in [`notebooks/06_calculate_afmultimer_predictions.ipynb`](./notebooks/06_calculate_afmultimer_predictions.ipynb) and run with the [LazyAF Google Colab notebook](https://colab.research.google.com/drive/1j7WJLcUHTR8BrjkWDaU549rFk6X5Zu18), but paste in the modified code in this notebook. The modified code adds lines to calculate the average pLDDT and include in the resulting CSV files for each comparison. Change the text in each cell to point to your data and rename the output CSV file. If you don't need the pLDDT information just use the default notebook code.
6. Save both the results and analysis folders from the AF-multimer predictions and the LazyAF analysis, which selects the highest ranking prediction JSON file and copies it into the analysis folder. Downloading the folder from Google Drive will automatically create a .zip archive and save this in your results folder for backup.
7. Analyze AF-multimer results compared to D-script predictions with [`notebooks/07_analyze-dscript-afmultimer-comparisons.ipynb`](./notebooks/07_analyze-dscript-afmultimer-comparisons.ipynb)
8. Code for analysis and plotting is in the [`notebooks/08_full-filtered-serine-protease-analysis.ipynb`](./notebooks/08_full-filtered-serine-protease-analysis.ipynb)
9. We also calculate additional AFmultimer metrics with the [`colabfold_analysis.py`](https://github.com/walterlab-HMS/AF2multimer-analysis/blob/main/colabfold_analysis.py) script which is part of the AF2multimer-analysis toolkit that underlies the [predictomes.org](https://predictomes.org/) webserver. To use the script:
   1.  Download the raw AF2multimer results from Google drive and unzip them. Google will sometimes split up large downloads into multiple zip files and not all results files that go together are in the same subdirectories. Make sure for this script to work all the `.json` and `.pae` files together as well as a `.done.txt` file with the name of the comparison
   2.  Clone the directory with `git clone https://github.com/walterlab-HMS/AF2multimer-analysis`
   3.  Run the script with `python3 colabfold_analysis.py <AF2MULTIMER_RESULTS_DIRECTORY>`. This will create a new folder with the same name as your input directory with `_analysis` appended at the end. This creates three CSV files: a contacts, interfaces, and summary file. The `summary.csv` is the most useful for downstream steps and is how we calculated the pDockQ score for each interaction prediction. The contacts and interfaces files give more information for each residue the calculated metrics, whereas the summary file is an average across all residues.
10. Obtain the predicted expression of each human serine protease hit from the [NCBI Gene Database](https://www.ncbi.nlm.nih.gov/gene). For each gene, we downloaded primary tissue expression data generated by [Fagerberg et al. 2014](https://pubmed.ncbi.nlm.nih.gov/24309898/), which analyzed gene expression in 27 different tissues from 95 healthy humans. Code for analysis and plotting is found in [`notebooks/08_full-filtered-serine-protease-analysis.ipynb`](./notebooks/08_full-filtered-serine-protease-analysis.ipynb)

### Compute Specifications

Scripts and analysis run on MacOS except for AF-multimer runs which were through Google Colab. Most of the runs on Google Colab were with A100 GPUs, and when not available with either V100 or T4 GPUs.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide--credit-for-contributions.md)
