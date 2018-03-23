#!/usr/bin/env python
# Copyright (C) 2018 Sur Herrera Paredes

import fyrd
import argparse
import os


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Uses EMBOSS' transeq to translate all fasta "
                          "sequences in all fasta files in a directory.\n"
                          "Currently only handles standard genetic code.")

    # Define required arguments
    required.add_argument("--indir", help=("Directory where fasta files are"),
                          required=True, type=str)
    required.add_argument("--outdir", help=("Directory where to write output "
                                            "files"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--fasta_suffix", help=("Suffix of fasta files in "
                                                "indir"),
                        type=str, default='.fna')
    parser.add_argument("--out_suffix", help=("Suffix of output files"),
                        type=str, default='.faa')
    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()

    fasta_files = os.listdir(args.indir)
    fasta_files = list(filter(lambda f: f.endswith(args.fasta_suffix),
                              fasta_files))
    print(fasta_files)
