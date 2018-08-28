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
    parser.description = ("It creates a mapping file that groups genomes "
                          "into batches. This file can be used by Nextflow "
                          "to decide how to split the genomes for "
                          "processing.\n"
                          "It assumes that the directory passed has a set "
                          "of subdirectories corresponding to some arbitrary "
                          "group. Then each group subdirectory has one "
                          "subdirectory per genome where the only .fna "
                          "file in the genome subdirectory corresponds to "
                          "the contigs file")

    # Define required arguments
    required.add_argument("--indir", help=("Directory that containes genomes. "
                                           "See above for the expected "
                                           "directory structure."),
                          required=True, type=str)
    required.add_argument("--outfile", help=("Name of the mapping file to be "
                                             "created."),
                          required=True, type=str)
    required.add_argument("--batch_size", help=("Number of genomes per batch"),
                          required=True, type=int)

    # Define other arguments
    parser.add_argument("--outdir", help=("Directory where to create "
                                          "symbolic link structure"),
                        type=str, default='genome_links/')

    parser.add_argument("--create_symlinks", help=("Flag indicating whether "
                                                   "to create a symbolic link "
                                                   "directory structure "
                                                   "for the defined batches"),
                        type=str,
                        default=False,
                        action='store_true')

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = process_arguments()
