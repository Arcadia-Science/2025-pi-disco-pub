#! /usr/bin/env python3

import argparse

from Bio import SeqIO


def filter_fasta_by_locus_tags(fasta_file, locus_tags_file, output_file):
    # Read the locus tags into a set for quick lookup
    with open(locus_tags_file) as lt_file:
        locus_tags = set(lt_file.read().splitlines())

    # Filter the FASTA records and write matching ones to the output file
    with open(output_file, "w") as out_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Check if any of the locus tags is a substring of the record ID
            if any(locus_tag in record.id for locus_tag in locus_tags):
                SeqIO.write(record, out_file, "fasta")


def main():
    parser = argparse.ArgumentParser(description="Filter FASTA file by locus tags.")
    parser.add_argument("fasta_file", type=str, help="Path to the input FASTA file.")
    parser.add_argument("locus_tags_file", type=str, help="Path to locus tags file.")
    parser.add_argument("output_file", type=str, help="Path for the output FASTA file.")

    args = parser.parse_args()

    filter_fasta_by_locus_tags(args.fasta_file, args.locus_tags_file, args.output_file)


if __name__ == "__main__":
    main()
