#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

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
import pandas

def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Script that takes a list of PATRIC genome "
                          "IDs and downloads the genomes to a specified "
                          "location.\n"
                          "It can also take a grouping factor per genome "
                          "and organize the downloads in subdirectories "
                          "according to that grouping factor. Spaces in "
                          "the name will be converted to underscores.\n"
                          "A name can also be passed as a column. This "
                          "name will be used to name the subdirectory "
                          "containing each genome data. Spaces in the name "
                          "will be converted to underscores.")

    # Define required arguments
    required.add_argument("--genomes", help=("Text file specifying the "
                                             "genome IDs to download. One "
                                             "genome per line only. It can "
                                             "have multiple tab-delimited "
                                             "columns in which case the "
                                             "columns corresponding to the "
                                             "and group need to be specified"),
                          required=True, type=str)
    required.add_argument("--outdir", help=("Directory where to download the "
                                            "genomes. The script will create "
                                            "one subdirectory per genome. A "
                                            "name for the subdirectory can "
                                            "also be specified in a column "
                                            "of the genomes file."))

    # Define other arguments
    parser.add_argument("--id_col", help=("Column number where the genome "
                                          "ID is stored. If not passed, "
                                          "the script will assume that the "
                                          "whole line is the genome ID."),
                        default=0,
                        type=int)
    parser.add_argument("--group_col", help=("Column number where the genome "
                                             "group is stored. If passed, "
                                             "id_col must also be specified. "
                                             "If passed, the script will "
                                             "create a subdirectory per group "
                                             "and download each genome into "
                                             "the corresponding "
                                             "subdirectory."),
                        default=0,
                        type=int)
    parser.add_argument("--name_col", help=("Column number where the genome "
                                            "name is stored. If passed, "
                                            "id_col must also be specified. "
                                            "If passed, the script will "
                                            "name the subdirectory per genome "
                                            "according to the specified "
                                            "column; otherwise, the genome ID "
                                            "will be used to name the "
                                            "subdirectory"),
                        default=0,
                        type=int)
    parser.add_argument("--header", help=("Flag indicating if genomes file "
                                          "contains headers."),
                        default=False,
                        action="store_true")
    parser.add_argument("--overwrite", help=("Flag indicating if "
                                             "subdirectories for genomes that "
                                             "already exist must be "
                                             "overwritten."),
                        default=False,
                        action="store_true")
    parser.add_argument("--failed", help=("Name of the file where to write "
                                          "failed downloads."),
                        default="failed.txt",
                        type=str)

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    if args.name_col > 0 and args.id_col <= 0:
        raise argparse.ArgumentError("--id_col is required if --name_col "
                                     "is passed")

    if args.group_col > 0 and args.id_col <= 0:
        raise argparse.ArgumentError("--id_col is required if --group_col "
                                     "is passed")

    return args
