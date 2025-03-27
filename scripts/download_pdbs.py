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
    upload_results_to_s3,
    read_sequences_from_file,
    replace_x_with_alanine,
)

import aiohttp
import argparse
import asyncio
import os
import sys

rate_limiter = RateLimiter(max_calls=100, period=1)  # 100 calls per second


async def is_job_finished(session: aiohttp.ClientSession, job_id: str, seq_id: str):
    """
    Checks if a job is finished asyncronously.

    Args:
        session (aiohttp.ClientSession): The aiohttp ClientSession.
        job_id (str): The id of the job to check.
        seq_id (str): The id (name) of the sequence for which the job was submitted.

    Returns:
        bool: True if the job is finished, False otherwise. If True, it'll save the file
    """

    for i in range(MAX_RETRIES):
        try:
            with rate_limiter:
                async with session.get(
                    f"{API_URL}/result/{job_id}", headers=HEADERS
                ) as response:
                    response.raise_for_status()
                    result = (await response.json())["result"]
                    if result != "not ready":
                        with open(
                            os.path.join(OUTPUT_DIR, f"{seq_id}.pdb"), "w"
                        ) as file:
                            file.write(result)
                        return True
                    else:
                        return False
        except aiohttp.ClientResponseError as e:
            if e.status in STATUS_FORCELIST:
                await asyncio.sleep(BACKOFF_FACTOR * (2**i))  # exponential backoff
            else:
                print(f"Error checking job {job_id}: {e}")
                return False
        except Exception as e:
            print(f"Error checking job {job_id}: {e}")
            return False
    print(
        f"Failed to check job {job_id} for sequence {seq_id} after {MAX_RETRIES} attempts."
    )
    return False


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
    successful_sequences_tsv: str, failed_sequences_tsv: str, total_sequences: int
) -> None:
    """
    Performs quality control checks on the processed sequences.

    Args:
        successful_sequences_tsv (str): The path to the TSV file containing successful downloads.
        failed_sequences_tsv (str): The path to the TSV file containing failed downloads.
        total_sequences (int): The total number of sequences processed.

    Returns:
        None
    """

    # Open the TSV files and count the lines (excluding headers)
    with open(successful_sequences_tsv, "r") as file:
        successful_lines = sum(1 for line in file) - 1
    with open(failed_sequences_tsv, "r") as file:
        failed_lines = sum(1 for line in file) - 1

    # Confirm the line numbers sum up to the total number of sequences
    assert (
        successful_lines + failed_lines == total_sequences
    ), "Line numbers do not match total sequences"

    assert failed_lines == 0, "There are failed downloads"

    # Check each PDB file in the successful_sequences_tsv has at least 10 rows of data
    with open(successful_sequences_tsv, "r") as file:
        next(file)  # Skip the header
        for line in file:
            pdb_name = line.split("\t")[2].strip()
            pdb_path = os.path.join(OUTPUT_DIR, pdb_name)
            with open(pdb_path, "r") as pdb_file:
                assert (
                    sum(1 for _ in pdb_file) >= 10
                ), f"PDB file {pdb_name} has less than 10 rows of data"


async def main(args=None):
    args = parse_args(args)

    output_dir = OUTPUT_DIR
    if not args.output_directory.startswith("s3://"):
        output_dir = args.output_directory

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create a sanitized job name
    pickle_file_name = get_metadata_file_name(args.job_name, "pkl", output_dir)

    # Load job_ids from checkpoint file if it exists
    job_ids = load_submitted_jobs(pickle_file_name)

    non_null_job_ids = {
        seq_id: job_id for seq_id, job_id in job_ids.items() if job_id is not None
    }

    to_download_job_ids = {}
    for seq_id, job_id in non_null_job_ids.items():
        if not os.path.exists(os.path.join(output_dir, f"{seq_id}.pdb")):
            to_download_job_ids[seq_id] = job_id

    async with aiohttp.ClientSession() as session:
        tasks = [
            is_job_finished(session, job_id, seq_id)
            for seq_id, job_id in to_download_job_ids.items()
        ]
        await asyncio.gather(*tasks)

    processed_sequences = {}
    failed_sequences = {}
    for seq_id, job_id in non_null_job_ids.items():
        pdb_file_path = os.path.join(output_dir, f"{seq_id}.pdb")
        if os.path.exists(pdb_file_path):
            processed_sequences[seq_id] = job_id
        else:
            failed_sequences[seq_id] = job_id

    fasta_file = download_fasta_file(args.fasta_file, output_dir)
    sequences = read_sequences_from_file(fasta_file)

    # Create a TSV file that includes the sequences that successfully downloaded
    successful_sequences_tsv = get_metadata_file_name(
        f"{args.job_name}_succesful_download", "tsv", output_dir
    )
    with open(successful_sequences_tsv, "w") as file:
        file.write(
            "gene_name\tfasta_header\tpdb_name\tis_original_sequence_valid\txtoa_replacement_count\tsequence\tsanitized_sequence\n"
        )
        for seq_id, sequence, header in sequences:
            if seq_id in processed_sequences:
                is_valid, sanitized_sequence = is_valid_amino_acid_sequence(sequence)
                if is_valid:
                    file.write(
                        f"{seq_id}\t{header}\t{seq_id}.pdb\t{is_valid}\t0\t{sequence}\t{sanitized_sequence}\n"
                    )
                else:
                    # Count the number of X's in the string
                    x_count = sanitized_sequence.count("X")

                    file.write(
                        f"{seq_id}\t{header}\t{seq_id}.pdb\t{is_valid}\t{x_count}\t{sequence}\t{replace_x_with_alanine(sanitized_sequence)}\n"
                    )

    # Create a TSV file that includes the sequences that failed to download
    failed_sequences_tsv = get_metadata_file_name(
        f"{args.job_name}_failed_download", "tsv", output_dir
    )
    with open(failed_sequences_tsv, "w") as file:
        file.write("gene_name\tjob_id\n")
        for seq_id, job_id in failed_sequences.items():
            file.write(f"{seq_id}\t{job_id}\n")

    qc_check(successful_sequences_tsv, failed_sequences_tsv, len(non_null_job_ids))

    # Upload results and the TSV to S3, if applicable
    upload_results_to_s3(args.output_directory, output_dir)


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
