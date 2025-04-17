### ESMFold API scripts

These are a couple utility scripts that submit folding jobs to the Arcadia Science hosted ESMFold API and to download the folded proteins. `enqueue.py` takes the FASTA file output (provided from S3), folds them using the Arcadia Science hosted ESMFold API. Example usage:

```
python scripts/enqueue.py --fasta-file s3://arcadia-pi-disco/tick_PUFs_PIs_1200.fasta --job-name "09/12 PI discovery" --output-directory s3://arcadia-pi-disco/pdb-files
```

As a side-effect, it creates a pickle file that stores the `sequence_id` (from the FASTA file) to the `job_id` (job ID used to fetch the folded results from the ESMFold API) and a TSV file that includes the gene names for failed job submissions (ends with `_submission_errors.tsv`).

Inputs:

- `--fasta-file` specifies the input FASTA file to be used.
- `--job-name` specifies the prefix to be used for the metadata files
- `--output-directory` species the local or remote (S3 URIs work) directory to put the results in. If an S3 URI is specified, the local results will be stored in the `outputs/` directory.

The second script `download_pdbs.py` uses the pickle from the other script to start downloading PDB files specified by the `--output-directory`. Example usage:

```
python scripts/download_pdbs.py --fasta-file s3://arcadia-pi-disco/tick_PUFs_PIs_1200.fasta --job-name "09/12 PI discovery" --output-directory s3://arcadia-pi-disco/pdb-files
```

As a side-effect, it creates two TSV files: one with all the successfully downloaded sequences and one with all the failed to download ones

Inputs:

- `--fasta-file` specifies the input FASTA file to be used.
- `--job-name` specifies the prefix to be used for the metadata files
- `--output-directory` species the local or remote (S3 URIs work) directory to put the results in. If an S3 URI is specified, the local results will be stored in the `outputs/` directory.
