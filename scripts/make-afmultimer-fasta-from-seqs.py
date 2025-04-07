import argparse
import os
from itertools import product

from Bio import SeqIO


def read_fasta(file_path):
    """
    Reads a FASTA file and returns a dictionary of sequences with the header as the key.
    """
    with open(file_path) as file:
        sequences = {record.id: str(record.seq) for record in SeqIO.parse(file, "fasta")}
    return sequences


def get_uniprot_accession(candidate_id):
    """
    Extracts the UniProt accession number from a candidate identifier.
    Assumes the format is like 'sp|A0A1B0GVH4|PRS51_HUMAN'.
    """
    parts = candidate_id.split("|")
    if len(parts) == 3:
        return parts[1]  # The UniProt accession number is the second element.
    return candidate_id  # Return the original ID if the expected format is not found.


def generate_combinations(baits, candidates, output_dir, length_threshold):
    """
    Generates all combinations of bait and candidate sequences,
    calculates the length of the combined sequence,
    and writes each combination to a separate FASTA file in a subdirectory based on sequence length.
    """
    under_threshold_dir = os.path.join(output_dir, f"under_{length_threshold}")
    over_threshold_dir = os.path.join(output_dir, f"over_{length_threshold}")
    # Ensure subdirectories exist
    os.makedirs(under_threshold_dir, exist_ok=True)
    os.makedirs(over_threshold_dir, exist_ok=True)

    for bait_id, candidate_id in product(baits, candidates):
        uniprot_accession = get_uniprot_accession(candidate_id)
        formatted_header = f"[{bait_id}_vs_{uniprot_accession}]"
        combined_sequence = baits[bait_id] + ":" + candidates[candidate_id]

        # Choose the directory based on sequence length
        if len(combined_sequence) <= length_threshold:
            output_file = os.path.join(
                under_threshold_dir, f"{bait_id}_vs_{uniprot_accession}.fasta"
            )
        else:
            output_file = os.path.join(
                over_threshold_dir, f"{bait_id}_vs_{uniprot_accession}.fasta"
            )

        with open(output_file, "w") as file:
            file.write(f">{formatted_header}\n")
            file.write(f"{combined_sequence}\n")


def main(bait_file, candidate_file, output_dir, length_threshold):
    baits = read_fasta(bait_file)
    candidates = read_fasta(candidate_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    generate_combinations(baits, candidates, output_dir, length_threshold)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate combinations of bait/candidate proteins and sort by sequence length."
    )
    parser.add_argument("bait_file", help="FASTA file containing bait proteins.")
    parser.add_argument("candidate_file", help="FASTA file containing candidate proteins.")
    parser.add_argument("output_dir", help="Directory to store output FASTA files.")
    parser.add_argument(
        "length_threshold", type=int, help="Length threshold to separate the directories."
    )

    args = parser.parse_args()

    main(args.bait_file, args.candidate_file, args.output_dir, args.length_threshold)
