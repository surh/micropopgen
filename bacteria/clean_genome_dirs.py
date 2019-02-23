#!/usr/bin/env python

# Copyright (C) 2019 Sur Herrera Paredes

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
import shutil


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Script that runs throuhg a directory structure "
                          "containing genomes, and creates a parallel "
                          "structure in a specified folder, discarding "
                          "unwanted genomes")

    # Define required arguments
    required.add_argument("indir", help=("Directory with genome direcetories"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--keep", help=("File with genomes to keep, one per "
                                        "line"),
                        type=str)
    parser.add_argument("--discard", help=("File with genomes to discard, "
                                           "one per line"),
                        type=str)
    parser.add_argument("--outdir", help=("Directory where to place cleaned "
                                          "genomes."),
                        type=str,
                        default="output/")
    parser.add_argument("--copy", help=("Copy wanted genomes instead of link"),
                        default=False, action="store_true")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    if args.keep is None and args.discard is None:
        raise ValueError("At least one of keep and discard must be provided")

    return args


def list_level_dirs(path, level=1):
    """Recursive function to return all directories at a specified level"""

    if level == 1:
        # Base case, when current level is wanted
        dirs = os.listdir(path)
    elif level > 1:
        # If we need to go deeper
        next_level_dirs = os.listdir(path)
        dirs = []
        for d in next_level_dirs:
            tempd = list_level_dirs(os.path.join(path, d), level - 1)
            dirs.extend([os.path.join(d, t) for t in tempd])
    else:
        # Error
        raise ValueError("ERROR: level must be a positive integer")

    return dirs


def read_list(path):
    """Read file into list, one value per line"""

    with open(path, 'r') as ih:
        res = []
        for line in ih:
            res.append(line.rstrip("\n"))
    ih.close()

    return(res)


if __name__ == "__main__":
    args = process_arguments()

    keep = []
    discard = []
    if args.keep is not None:
        keep = read_list(args.keep)
    if args.discard is not None:
        discard = read_list(args.discard)
