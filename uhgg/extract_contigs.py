#!/usr/bin/env python
# Copyright (C) 2020 Sur Herrera Paredes

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from Bio import SeqIO
import os
import argparse


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Read a genome fasta file from UHGG and separate "
                          "into one record per contig.")

    # Define required arguments
    required.add_argument("--input",
                          help=("Path to input fasta file"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--outdir",
                        help=("Path to output directory to create"),
                        type=str,
                        default="output/")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


def split_fasta_records(input, outdir):
    """Split a fasta file into a new fasta file per record.
    Save new files in a previously exisiting directory"""

    i = 0
    for record in SeqIO.parse(input, 'fasta'):
        record_fasta = outdir + '/' + record.id + '.fasta'
        SeqIO.write(record, record_fasta, 'fasta')
        i = i + 1

    return i


if __name__ == "__main__":
    args = process_arguments()

    if os.path.isdir(args.outdir):
        raise FileExistsError("Output directory already exists")
    else:
        os.mkdir(args.outdir)

    n_recs = split_fasta_records(args.input, args.outdir)
    print("{} records found.".format(n_recs))
