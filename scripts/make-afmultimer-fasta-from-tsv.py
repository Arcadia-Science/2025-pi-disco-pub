#! /usr/bin/env python3

import argparse
import csv
import os

from Bio import SeqIO


def extract_protein_id(candidate):
    """Extract the protein ID from the candidate column,
    handling both [protein_id=XXX] and | formats."""
    if "[protein_id=" in candidate:
        return candidate.strip("[]").split("=")[1]
    else:
        parts = candidate.split("|")
        # Assuming the relevant ID is always after the first |, typical in UniProt formats
        return parts[1] if len(parts) > 1 else candidate


def parse_tsv(tsv_file_path):
    """Parse the TSV file to map bait proteins to their
    candidate entries, adapting to different candidate formats."""
    entries_to_combine = {}
    with open(tsv_file_path) as tsv_file:
        for row in csv.reader(tsv_file, delimiter="\t"):
            bait, candidate_with_prefix = row[0], row[1]
            candidate = extract_protein_id(candidate_with_prefix)
            entries_to_combine.setdefault(bait, set()).add(candidate)
    return entries_to_combine


def parse_fasta(fasta_file_path):
    """Parse the FASTA file and return a dictionary of sequences indexed by their headers."""
    fasta_sequences = {}
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        fasta_sequences[record.id] = record
    return fasta_sequences


def find_sequence_by_id(fasta_sequences, search_id):
    """Find a sequence in the fasta_sequences dictionary by partial ID match."""
    for seq_id in fasta_sequences:
        if search_id in seq_id:
            return fasta_sequences[seq_id].seq
    return None


def create_fasta_files(entries_to_combine, fasta_sequences, output_directory):
    """Create new FASTA files combining bait and candidate sequences."""
    processed_count = 0
    for bait, candidates in entries_to_combine.items():
        bait_sequence = find_sequence_by_id(fasta_sequences, bait)
        if not bait_sequence:
            print(f"Warning: No match found for bait {bait}")
            continue

        for candidate in candidates:
            candidate_sequence = find_sequence_by_id(fasta_sequences, candidate)
            if not candidate_sequence:
                print(f"Warning: No match found for candidate {candidate}")
                continue

            combined_sequence = str(bait_sequence) + ":" + str(candidate_sequence)
            output_file_path = os.path.join(output_directory, f"{bait}_vs_{candidate}.fasta")
            with open(output_file_path, "w") as output_file:
                output_file.write(f">[{bait}_{candidate}]\n{combined_sequence}\n")
            processed_count += 1

    return processed_count


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine FASTA sequences based on TSV file entries"
    )
    parser.add_argument("tsv_file_path", help="Path to the TSV file")
    parser.add_argument("fasta_file_path", help="Path to the FASTA file")
    parser.add_argument("output_directory", help="Directory to save the output FASTA files")

    args = parser.parse_args()

    # Check if output directory exists; if not, create it
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    entries_to_combine = parse_tsv(args.tsv_file_path)
    fasta_sequences = parse_fasta(args.fasta_file_path)
    processed_count = create_fasta_files(entries_to_combine, fasta_sequences, args.output_directory)

    print(f"{processed_count} FASTA files have been created.")
