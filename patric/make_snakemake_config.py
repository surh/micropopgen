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

import argparse


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Creates a config file that contains the "
                          "names of the expected output files.")

    # Define required arguments
    required.add_argument("--infile", help=("A tab delimited file "
                                            "indicating species and genome "
                                            "id."),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--outfile", help=("Name of the output file"),
                        type=str,
                        default="config.yml")
    parser.add_argument("--species_col", help=("Columb wheere species is "
                                               "stored"),
                        type=int, default=1)
    parser.add_argument("--id_col", help=("Column where genome idis stored."),
                        type=int, default=2)
    parser.add_argument("--header", help=("Flag indicating if header is "
                                          "present"),
                        default=False, action="store_true")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    return args


def make_snakemake_config(infile, outfile,
                          species_col=1,
                          id_col=2,
                          header=False):
    """Read tab delimited file and create .yaml file which
    has the list of expected outputs from the snakemake pipeline"""

    id_col = id_col - 1
    species_col = species_col - 1
    with open(infile, 'r') as ih, open(outfile, 'w') as oh:
        oh.write("files:\n")
        print("Reading list of species")
        i = 0
        for line in ih:
            line = line.rstrip()
            if i == 0 and header:
                i = i + 1
                continue

            LINE = line.split("\t")
            species = LINE[species_col]
            id = LINE[id_col]

            filename = species + '-' + id + '.gff'
            oh.write("\t" + id + ': ' + filename + "\n")
            i = i + 1
    print("Done")
    ih.close()
    oh.close()

    return i


if __name__ == "__main__":
    args = process_arguments()

    make_snakemake_config(infile=args.infile,
                          outfile=args.outfile,
                          species_col=args.species_col,
                          id_col=args.id_col,
                          header=args.header)
