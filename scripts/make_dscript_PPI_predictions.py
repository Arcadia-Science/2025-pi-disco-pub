import argparse
import subprocess
from itertools import product

import pandas as pd


def read_fasta(fasta_file):
    with open(fasta_file) as file:
        lines = file.readlines()
    headers = []
    sequences = []
    seq = ""
    for line in lines:
        if line.startswith(">"):
            if seq:
                sequences.append(seq)
                seq = ""
            # Remove the '>' character and any leading/trailing whitespace
            headers.append(line.strip().lstrip(">"))
        else:
            seq += line.strip()
    if seq:
        sequences.append(seq)
    return headers, sequences


def write_tsv(bait_headers, target_headers, filename):
    pairs = list(product(bait_headers, target_headers))
    df = pd.DataFrame(pairs)
    df.to_csv(filename, sep="\t", index=False, header=False)


def combine_fastas(fasta1, fasta2, combined_filename):
    with open(fasta1) as file1, open(fasta2) as file2, open(combined_filename, "w") as outfile:
        for line in file1:
            outfile.write(line)
        for line in file2:
            outfile.write(line)


def call_dscript(tsv_file, fasta_file, model_path):
    command = f"dscript predict --pairs {tsv_file} --seqs {fasta_file} --model {model_path}"
    process = subprocess.run(command, shell=True, text=True)
    if process.stderr:
        print("Error:", process.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Create TSV of bait-prey pairs, combine FASTA files, and run dscript."
    )
    parser.add_argument("bait_fasta", help="Path to the bait FASTA file")
    parser.add_argument("target_fasta", help="Path to the FASTA file of target candidates")
    parser.add_argument("output_tsv", help="Path to the output TSV file")
    parser.add_argument("combined_fasta", help="Path to the output combined FASTA file")
    parser.add_argument("model_path", help="Path to the dscript model file")

    args = parser.parse_args()

    bait_headers, _ = read_fasta(args.bait_fasta)
    target_headers, _ = read_fasta(args.target_fasta)

    # Create TSV of all combinations
    write_tsv(bait_headers, target_headers, args.output_tsv)

    # Combine FASTA files
    combine_fastas(args.bait_fasta, args.target_fasta, args.combined_fasta)

    # Call dscript
    call_dscript(args.output_tsv, args.combined_fasta, args.model_path)

    print("TSV and combined FASTA file created, and dscript predict command executed successfully.")


if __name__ == "__main__":
    main()
