from constants import (
    API_URL,
    BEARER_TOKEN,
    HEADERS,
    OUTPUT_DIR,
    BACKOFF_FACTOR,
    MAX_RETRIES,
    STATUS_FORCELIST,
)
from ratelimiter import RateLimiter
from typing import Optional, Sequence
from utils import (
    download_fasta_file,
    get_metadata_file_name,
    is_valid_amino_acid_sequence,
    load_submitted_jobs,
    read_sequences_from_file,
    replace_x_with_alanine,
)

import aiohttp
import argparse
import asyncio
import os
import pickle
import sys

"""
Usage:
mamba env create --file env.yml
mamba activate pi-disco

python 03_run_esmfold/enqueue.py --fasta-file s3://arcadia-pi-disco/tick_PUFs_PIs_1200.fasta --job-name "09/12 PI discovery" --output-directory s3://arcadia-pi-disco/pdb-files
"""

rate_limiter = RateLimiter(max_calls=100, period=1)  # 100 calls per second


async def submit_job(session: aiohttp.ClientSession, sequence: str):
    """
    Submits a job to the API asynchronously.

    Args:
        session (aiohttp.ClientSession): The aiohttp ClientSession.
        sequence (str): The sequence to be submitted.

    Returns:
        str: The job id.
    """

    for i in range(MAX_RETRIES):
        try:
            with rate_limiter:
                async with session.post(
                    f"{API_URL}/fold",
                    headers=HEADERS,
                    json={"protein_sequence": sequence},
                ) as response:
                    response.raise_for_status()
                    return (await response.json())["job_id"]
        except aiohttp.ClientResponseError as e:
            if e.status in STATUS_FORCELIST:
                await asyncio.sleep(BACKOFF_FACTOR * (2**i))  # exponential backoff
            else:
                print(f"Error submitting job for sequence {sequence}: {e}")
                return None
        except Exception as e:
            print(f"Error submitting job for sequence {sequence}: {e}")
            return None
    print(f"Failed to submit job for sequence {sequence} after {MAX_RETRIES} attempts.")
    return None


def parse_args(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """
    Parses command line arguments.

    Args:
        args (Optional[Sequence[str]]): The command line arguments. If None, the arguments are taken from sys.argv.

    Returns:
        argparse.Namespace: The parsed arguments.
    """
    Description: str = "Submit protein folding jobs."
    Epilog: str = "Example usage: main.py <FASTA_FILE> <JOB_NAME>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-f", "--fasta-file", type=str, required=True, help="Path to the FASTA file."
    )
    parser.add_argument(
        "-j", "--job-name", type=str, required=True, help="Name for the folding job."
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        type=str,
        required=True,
        help="Name for the output directory (could be a local path or an S3 URI).",
    )
    return parser.parse_args(args)


def qc_check(
    failed_sequences_tsv: str, successful_sequences: int, total_sequences: int
) -> None:
    """
    Performs quality control checks on the processed sequences.

    Args:
        failed_sequences_tsv (str): The path to the TSV file containing failed sequences.
        successful_sequences (int): The total number of successfully submitted sequences
        total_sequences (int): The total number of sequences processed.

    Returns:
        None
    """
    # Open the TSV files and count the lines (excluding headers)
    with open(failed_sequences_tsv, "r") as file:
        failed_lines = sum(1 for line in file) - 1

    assert (
        successful_sequences + failed_lines == total_sequences
    ), "Line numbers do not match total sequences"
    assert failed_lines == 0, "There are failed submissions"


async def main(args=None):
    args = parse_args(args)

    output_dir = OUTPUT_DIR
    if not args.output_directory.startswith("s3://"):
        output_dir = args.output_directory

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fasta_file = download_fasta_file(args.fasta_file, output_dir)
    sequences = read_sequences_from_file(fasta_file)
    total_sequence_count = len(sequences)

    # Create a sanitized job name
    pickle_file_name = get_metadata_file_name(args.job_name, "pkl", output_dir)

    # Load job_ids from checkpoint file if it exists
    job_ids = load_submitted_jobs(pickle_file_name)

    # Filter sequences based on is_valid_amino_acid
    valid_sequences = []
    for seq_id, sequence, header in sequences:
        is_valid, sanitized_sequence = is_valid_amino_acid_sequence(sequence)
        # While processing these results, we came across a couple amino acid
        # sequences that include the X single character code.
        # "X" is not a valid AA - it is what happens when there is a sequencing ambiguity
        #  in the underlying nucleotide (DNA or RNA sequence). These ambiguities will
        # be assigned an "N" instead of a real nucleotide (ATCG). An "N" will turn into a "X"
        # in the protein because there is no way to confidently turn it to an AA.
        # For the purposes of this experiment, we'll turn all "X" AAs into "A"s
        if is_valid:
            valid_sequences.append((seq_id, sanitized_sequence, header))
        else:
            valid_sequences.append(
                (seq_id, replace_x_with_alanine(sanitized_sequence), header)
            )

    # Only submit sequences if job_ids[seq_id] is False-y
    valid_sequences = [
        (seq_id, sequence, header)
        for seq_id, sequence, header in valid_sequences
        if not job_ids.get(seq_id)
    ]

    async with aiohttp.ClientSession() as session:
        tasks = [
            submit_job(session, sequence)
            for seq_id, sequence, header in valid_sequences
        ]
        job_ids_results = await asyncio.gather(*tasks)

        for (sequence, seq_id), job_id in zip(valid_sequences, job_ids_results):
            job_ids[seq_id] = job_id

        # Make sure to record the job_ids
        with open(pickle_file_name, "wb") as file:
            pickle.dump(job_ids, file)

    null_job_ids = {
        seq_id: job_id for seq_id, job_id in job_ids.items() if job_id is None
    }

    non_null_job_ids = {
        seq_id: job_id for seq_id, job_id in job_ids.items() if job_id is not None
    }

    # Create a TSV file that includes the sequences that failed to start
    failed_sequences_tsv = get_metadata_file_name(
        f"{args.job_name}_submission_errors", "tsv", output_dir
    )
    with open(failed_sequences_tsv, "w") as file:
        file.write("gene_name\n")
        for seq_id in null_job_ids.keys():
            file.write(f"{seq_id}\n")

    qc_check(failed_sequences_tsv, len(non_null_job_ids), total_sequence_count)


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
