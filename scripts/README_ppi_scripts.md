# Documentation for running command-line scripts

There are two scripts to make the input FASTA files for Alphafold-multimer. One takes in two FASTA files, a bait FASTA and candidate FASTA, and the script makes pairwise combinations of the bait:candidate sequences separated by `:`. The other takes in a TSV of specific combinations that are in a single input FASTA file.

Additionally there is a length threshold argument in the first script to separate the output FASTA sequences by sequences longer than the input threshold, for example something like 700 amino acids. This is if you want to only run AF-multimer predictions for sequences underneath a certain length if you are trying to avoid long compute times/limited compute resources. The length of a combined sequence complex and the amount of compute time can depend on (anecdotally) the quality of the structures plus the type of GPUs available (A100 are about 10X faster than T4 GPUs but are more expensive). So depending on this information is how you should choose the length threshold at which to split sets of proteins.

```
usage: make-afmultimer-fasta-from-seqs.py [-h]
                                        bait_file candidate_file output_dir
                                        length_threshold
```

The script that takes in a TSV file of specific combinations to make

```
usage: make-afmultimer-fasta-from-tsv.py [-h]
                                         tsv_file_path fasta_file_path
                                         output_directory
```

The `make_dscript_PPI_predictions.py` script automates D-script predictions by inputting separate bait and candidate FASTA files. You will need to activate the conda environment with the pip installed D-script python package to use this script, documented in the main README.

```
usage: make_dscript_PPI_predictions.py [-h]
                                       bait_fasta target_fasta output_tsv
                                       combined_fasta model_path
```

R script is an analysis script for the D-script results.
