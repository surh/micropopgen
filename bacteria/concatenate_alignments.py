#!/usr/bin/env python
# Copyright (C) 2018 Sur Herrera Paredes
# This file is based on the AMPHORA
# pipeline by Martin Wu.

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
import os
import align_markers as am
from Bio import AlignIO


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Concatenate alns)

    # Define required arguments
    required.add_argument("--indir", help=("Input dir"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--output", help=("Output dir"),
                        type=str,
                        default="concatenated.aln.fasta")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()

    aln_files = os.listdir(args.indir)
    Alns = []
    for f in aln_files:
        a = AlignIO.read(os.path.join(args.indir, f), 'fasta')
        Alns.append(a)

    aln = am.concatenate_alignments(alns=Alns)
    AlignIO.write(aln, args.output, 'fasta')
