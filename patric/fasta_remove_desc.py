#!/usr/bin/env python
# Copyright (C) 2018 Sur Herrera Paredes

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

import argparse


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Script that takes a fasta file and removes "
                          "the description element from the header, "
                          "leaving only the IDs.")

    # Define required arguments
    required.add_argument("--infile", help=("Path to fasta input file"),
                          required=True, type=str)
    required.add_argument("--outfile", help=("Path to fasta output file"),
                          type=str,
                          required=True)

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    return args


def fasta_remove_desc(infile, outfile):
    """Take a fasta file and remove the description from
    every header"""

    with open(infile, 'r') as ih, open(outfile, 'w') as oh:
        print("Processing fasta file")
        for line in ih:
            line = line.rstrip()
            if line.startswith('>'):
                id = line.split()[0]
                oh.write(id + '\n')
            else:
                oh.write(line + '\n')
    print("Done")

    return


if __name__ == "__main__":
    args = process_arguments()

    fasta_remove_desc(args.infile, args.outfile)
