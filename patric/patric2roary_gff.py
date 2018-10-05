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
    parser.description = ("Takes a PATRIC GFF 3 file and changes edits it "
                          "to simplofy the output")

    # Define required arguments
    required.add_argument("--infile", help=("Path to input GFF 3 file"),
                          required=True, type=str)
    required.add_argument("--outfile", help=("Path to output GFF 3 file"),
                          required=True, type=str)
    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    return args


def gff_patric2roary(infile, outfile):
    """Take a gff file from PATRIC and edit the chromosome and gene
    ID names for simplifying roary's output"""

    with open(infile, 'r') as ih, open(outfile, 'w') as oh:
        print("Processing GFF file")
        for line in ih:
            line = line.rstrip()
            if line.startswith('#'):
                oh.write(line + '\n')
            elif line == '':
                oh.write('\n')
            else:
                # Edit chromosome ID and gene name
                LINE = line.split("\t")
                LINE[0] = LINE[0].replace('accn|', '')
                gene_id = LINE[8].split(';')[0]
                # gene_id = gene_id.replace('ID=fig|', '')
                gene_id = gene_id.replace('fig|', '')
                gene_id = gene_id + ';'
                LINE[8] = gene_id
                newline = '\t'.join(LINE)
                oh.write(newline + "\n")
    print("Done")

    return


if __name__ == "__main__":
    args = process_arguments()

    gff_patric2roary(args.infile, args.outfile)
