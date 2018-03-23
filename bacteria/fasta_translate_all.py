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
    parser.add_argument("--transeq", help=("Binary executable of transeq "
                                           "tool of the EMBOSS package"),
                        type=str, default='transeq')
    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


def which(program):
    """Check if executable exists"""

    # From https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    # Under MIT license

    def is_exe(fpath):
        """Check if path is executable"""

        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


if __name__ == "__main__":
    args = process_arguments()

    # Check if exectuable exists
    if which(args.transeq) is None:
        raise FileNotFoundError("Executable for transeq not found")

    # Get list of fasta files from indir
    fasta_files = os.listdir(args.indir)
    fasta_files = list(filter(lambda f: f.endswith(args.fasta_suffix),
                              fasta_files))
    print(fasta_files)
