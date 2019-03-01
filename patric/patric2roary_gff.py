#!/usr/bin/env python
# Copyright (C) 2018-2019 Sur Herrera Paredes

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
from Bio import SeqIO
import re


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    # parser_format = argparse.RawTextHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Takes a PATRIC features.tab file and an .fna "
                          "file, and it creates a gff file that is compatible "
                          "with roary. "
                          "It makes no checks on the fna and features files.")

    # Define required arguments
    required.add_argument("--fna", help=("Path to input .fna file."),
                          required=True, type=str)
    required.add_argument("--features", help=("Path to input .features.tab "
                                              "file from PATRIC."),
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

    return

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


def patric_features_to_roary(fna_file, features_file, outfile):
    """Takes an fna file and a PATRIC features.tab file and creates
    a GFF3 file that is compatible with roary"""

    with open(features_file, 'r') as feat, open(outfile, 'w') as out:
        print("Processsing features file")
        feat.readline()
        oh.write("##gff-version 3\n")
        for line in feat:
            # Read line
            line = line.rstrip()
            LINE = line.split("\t")
            seqid = LINE[2]
            source = LINE[3]
            type = LINE[4]
            start = LINE[9]
            end = LINE[10]
            strand = LINE[11]
            phase = '0'
            id = re.sub('\w+\|', '', LINE[5])

            # Create attributes
            attr = "ID={0};locus_tag={0}".format(id)

            # Create new line
            gff_line = '\t'.join([seqid, source, type, start, end, strand,
                                  phase, attr])
            out.write(''.join([gff_line, '\n']))
        feat.close()

    return


if __name__ == "__main__":
    args = process_arguments()

    # gff_patric2roary(args.infile, args.outfile)
    patric_features_to_roary(args.fna, args.features, args.outfile)
