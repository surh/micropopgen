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


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Script that takes a directory of marker files "
                          "and makes an alignment per marker file and "
                          "filters the alignment. Based on the AMPHORA "
                          "pipeline")

    # Define required arguments
    required.add_argument("--indir", help=("Directory with the individual "
                                           "marker files"),
                          required=True, type=str)
    required.add_argument("--outdir", help=("Directory where to store the "
                                            "output"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--muscle", help=("Executable of muscle"),
                        type=str,
                        default="muscle")
    parser.add_argument("--marker_suffix", help=("Suffix of files with "
                                                 "sequences"),
                        type=str, default='.faa')

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()

    # Check input files
    if not os.path.isdir(args.indir):
        raise FileNotFoundError("Indir not found")

    # Create output directory
    if os.path.isdir(args.outdir):
        raise FileExistsError("Outdir already exists")
    else:
        os.mkdir(args.outdir)
