from Bio import SeqIO
from constants import AWS_ACCESS_KEY_ID, AWS_REGION, AWS_SECRET_ACCESS_KEY

import boto3
import os
import pickle
import re


def sanitize_name(name: str) -> str:
    """
    Sanitizes a name by removing non-alphanumeric characters and converting to lowercase.

    Args:
        name (str): The name to sanitize.

    Returns:
        str: The sanitized and lowercased name.
    """
    name = name.replace(" ", "_")
    return re.sub(r"\W+", "", name).lower()


def get_metadata_file_name(job_name: str, extension: str, output_dir: str) -> str:
    """
    Generates a metadata file name based on the job name.

    Args:
        job_name (str): The name of the job.
        extension (str): The extension for the filename
        output_dir (str): The output director for the file

    Returns:
        str: The generated file name.
    """
    sanitized_job_name = sanitize_name(job_name)
    return os.path.join(output_dir, f"{sanitized_job_name}.{extension}")


def load_submitted_jobs(pickle_file_name: str) -> dict:
    """
    Loads submitted jobs from a pickle file.

    Args:
        pickle_file_name (str): The name of the pickle file.

    Returns:
        dict: A dictionary of job ids. If the file does not exist, an empty dictionary is returned.
    """
    try:
        with open(pickle_file_name, "rb") as file:
            job_ids = pickle.load(file)
    except FileNotFoundError:
        job_ids = {}
    return job_ids


def read_sequences_from_file(file_path: str):
    """
    Reads sequences from a file.

    Args:
        file_path (str): The path to the file.

    Returns:
        list: A list of tuples containing the sequence and its id.
    """
    sequences = [
        (record.id, str(record.seq), str(record.description))
        for record in SeqIO.parse(file_path, "fasta")
    ]
    return sequences


def is_valid_amino_acid_sequence(sequence: str):
    """
    Reads sequences from a file.

    Args:
        sequence (str): FASTA sequence for the protein

    Returns:
        tuple: A tuple of valid flag (True or False) and the sanitized sequence
    """

    if len(sequence) == 0:
        return (False, sequence)

    # Define the standard amino acids
    amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

    # Sanitize the sequence by removing whitespaces and converting to uppercase
    sanitized_sequence = sequence.replace(" ", "").upper()

    # Check if each character in the sequence is a valid amino acid
    for char in sanitized_sequence:
        if char not in amino_acids:
            return (False, sanitized_sequence)
    return (True, sanitized_sequence)


def download_fasta_file(file_path: str, output_dir: str):
    """
    Downloads a FASTA file from S3 if needed

    Args:
        file_path (str): The path to the file.
        output_dir (str): The output director for the file

    Returns:
        str: The path to the staged file.
    """

    # Check if FASTA file is hosted on S3
    if file_path.startswith("s3://"):
        # Download FASTA file from S3 using boto3
        s3 = boto3.resource(
            "s3",
            aws_access_key_id=AWS_ACCESS_KEY_ID,
            aws_secret_access_key=AWS_SECRET_ACCESS_KEY,
            region_name=AWS_REGION,
        )
        bucket_name, key = file_path[5:].split("/", 1)
        download_path = os.path.join(output_dir, key)

        s3.Bucket(bucket_name).download_file(key, download_path)
        return download_path
    else:
        return os.path.abspath(file_path)


def upload_results_to_s3(output_directory: str, local_directory: str):
    """
    Uploads results to S3 if needed

    Args:
        output_directory (str): The path to the output directory.
        local_directory (str): Where the files are locally.

    Returns:
        None
    """

    if output_directory.startswith("s3://"):
        S3 = boto3.client(
            "s3",
            aws_access_key_id=AWS_ACCESS_KEY_ID,
            aws_secret_access_key=AWS_SECRET_ACCESS_KEY,
            region_name=AWS_REGION,
        )
        bucket_name, s3_directory = output_directory[5:].split("/", 1)
        for file_name in os.listdir(local_directory):
            S3.upload_file(
                os.path.join(local_directory, file_name),
                bucket_name,
                f"{s3_directory}/{file_name}",
            )


def replace_x_with_alanine(sequence: str):
    """
    Replaces X with alanine in the sequence

    Args:
        sequence (str): FASTA sequence for the protein

    Returns:
        str: Sanitized FASTA sequence for the protein
    """
    return sequence.replace("X", "A")
