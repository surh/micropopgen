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
import pandas as pd
import os


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Script that takes a list of PATRIC genome "
                          "IDs and downloads the genomes to a specified "
                          "location."
                          "It can also take a grouping factor per genome "
                          "and organize the downloads in subdirectories "
                          "according to that grouping factor. Spaces in "
                          "the name will be converted to underscores."
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
                                            "of the genomes file."),
                          required=True, type=str)

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


def read_genomes_file(file, id_col=0, group_col=0, name_col=0, header=False):
    """Uses pandas to read a file that contains genome IDs"""

    # Read table and get relevant columns
    if id_col <= 0:
        genomes = pd.read_table(file, sep="lnds555sgsg58s", engine='python',
                                dtype='str')
        genomes.columns = ['ID']
        genomes = genomes.assign(Name=genomes.ID)
    else:
        dat = pd.read_table(file, sep="\t")
        columns = [id_col - 1]
        colnames = ['ID']

        # Find columns
        if group_col > 0:
            columns.append(group_col - 1)
            colnames.append('Group')

        if name_col > 0:
            columns.append(name_col - 1)
        else:
            columns.append(id_col - 1)
        colnames.append('Name')

        # Get data
        genomes = dat.iloc[:, columns].astype('str').copy()
        genomes.columns = colnames

    # Replace whitespaces
    genomes.Name.replace(" ", "_", regex=True, inplace=True)
    if 'Group' in genomes.columns:
        genomes.Group.replace(" ", "_", regex=True, inplace=True)

    return genomes


def prepare_outdir(outdir, overwrite=False):
    """Check if output directory exists, and if not
    create it or raise error"""

    if os.path.isdir(outdir) and overwrite:
        print("\tOutput directory ({}) already exists".format(outdir))
    elif os.path.isdir(outdir) and not overwrite:
        raise FileExistsError("\tOutput directory {}"
                              " already exists".format(outdir))
    else:
        # If directory doesn't exist
        os.mkdir(outdir)

    return outdir


def download_genome_dir(id, name, outdir,
                        overwrite=False,
                        url="ftp://ftp.patricbrc.org/genomes/"):
    """Download genome directory from PATRIC database
    FTP server"""

    # Prepare output directory
    genome_dir = ''.join([outdir, '/', name])
    prepare_outdir(outdir=genome_dir, overwrite=overwrite)

    # Download to genome dir

    return


def download_genome_table(genomes, outdir, overwrite=False,
                          url="ftp://ftp.patricbrc.org/genomes/"):
    """Takes a pandas data frame where each row is a genomes
    and calls the function to download each one independently"""

    for i, r in genomes.iterrows():
        download_genome_dir(id=r['ID'], name=r['Name'],
                            outdir=outdir, overwrite=overwrite,
                            url="ftp://ftp.patricbrc.org/genomes/")

    return


if __name__ == "__main__":
    args = process_arguments()

    genomes = read_genomes_file(file=args.genomes,
                                id_col=args.id_col,
                                group_col=args.group_col,
                                name_col=args.name_col,
                                header=args.header)

    # Create main output Directory
    print("Creating overall output directory")
    prepare_outdir(outdir=args.outdir, overwrite=True)

    # Download genomes
    print("Downloadign genomes")
    if 'Group' in genomes:
        raise NotImplementedError("Group option is not implemented yet.")
    else:
        download_genome_table(genomes=genomes,
                              outdir=args.outdir,
                              overwrite=args.overwrite,
                              url="ftp://ftp.patricbrc.org/genomes/")
