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
import pandas as pd
import os
from ftplib import FTP
# import shutil


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
    parser.add_argument("--check", help=("Check whether an fna file was "
                                         "downloaded for each genome"),
                        default=False,
                        action="store_true")
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
        dat = pd.read_table(file, sep="\t", dtype={id_col - 1: str})
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
        print("\tCurrently at {}".format(os.getcwd()))
        print("\tCreating directory {}".format(outdir))
        os.mkdir(outdir)

    return outdir


def download_genome_dir(id, name, outdir,
                        overwrite=False,
                        url="ftp.patricbrc.org"):
    """Download genome directory from PATRIC database
    FTP server"""

    # Prepare output directory
    genome_dir = ''.join([outdir, '/', name])
    prepare_outdir(outdir=genome_dir, overwrite=overwrite)

    # Download to genome dir
    gdir = '/'.join(['genomes', id])
    cwd = os.getcwd()
    try:
        files = download_ftp_dir(ftp_url=url,  ftp_dir=gdir,
                                 ddir=genome_dir)
        success = True
    except:
        files = []
        success = False
        os.chdir(cwd)

    print("Downloaded {} files.".format(files))

    return success


def download_ftp_dir(ftp_url, ftp_dir, ddir):
    """Download a directoy via an FTP connection to an specified location.
    Rerurns list of files from ftp dir."""

    # Establish FTP connection and prepare download
    ftp = FTP(ftp_url)
    ftp.login()
    ftp.cwd(ftp_dir)
    ls = ftp.nlst()
    count = len(ls)
    print("found {} files".format(count))

    # Move to destination dir
    cwd = os.getcwd()
    os.chdir(ddir)

    # Download
    curr = 0
    for fn in ls:
        curr += 1
        print('Processing file {} ... {} of {} ...'.format(fn, curr, count))
        ftp.retrbinary('RETR ' + fn, open(fn, 'wb').write)

    # Close connection and go back
    ftp.quit()
    print("download complete.")
    os.chdir(cwd)

    return ls


def download_genome_table(genomes, outdir, overwrite=False,
                          url="ftp.patricbrc.org/genomes/"):
    """Takes a pandas data frame where each row is a genome
    and calls the function to download each one independently"""

    results = []
    for i, r in genomes.iterrows():
        print("====== Downloading genome {} ======".format(r['ID']))
        success = download_genome_dir(id=r['ID'], name=r['Name'],
                                      outdir=outdir, overwrite=overwrite,
                                      url="ftp.patricbrc.org")
        results.append([r['ID'], r['Name'], outdir, success])
        print("===================================")

    results = pd.DataFrame(results, columns=['ID', 'Name', 'Dir', 'Success'])

    return results


def check_genomes_dirs(indir):
    """Takes a directory that contains a number of genome subdirectories,
    and checks that every genome subdirtectory has a .fna file"""

    res = []
    if os.path.isdir(indir):
        # Find all genome subdirs
        specdirs = os.listdir(indir)
        for spec in specdirs:
            # Check if an fna file with the same name as the directory exists
            fna_filename = os.path.join(indir,
                                        spec,
                                        ''.join([spec, '.fna']))
            if os.path.isfile(fna_filename):
                success = True
            else:
                success = False

            res.append([spec, success])

        res = pd.DataFrame(res, columns=['ID', 'fna'])
    else:
        raise FileNotFoundError("Genome directory does not exist")

    return res


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
    print("Downloading genomes")
    if 'Group' in genomes:
        # raise NotImplementedError("Group option is not implemented yet.")
        results = pd.DataFrame()
        fna = pd.DataFrame()
        for n, g in genomes.groupby('Group'):
            print("------------------------------------------------")
            print("Processing group {}".format(n))
            group_dir = ''.join([args.outdir, '/', n])
            prepare_outdir(outdir=group_dir, overwrite=True)
            r = download_genome_table(genomes=g,
                                      outdir=group_dir,
                                      overwrite=args.overwrite,
                                      url="ftp.patricbrc.org")
            if(args.check):
                fna = fna.append(check_genomes_dirs(group_dir),
                                 ignore_index=True,
                                 sort=False)
            results = results.append(r)
            print("------------------------------------------------")
    else:
        results = download_genome_table(genomes=genomes,
                                        outdir=args.outdir,
                                        overwrite=args.overwrite,
                                        url="ftp.patricbrc.org")
        if(args.check):
            fna = fna.append(check_genomes_dirs(args.outdir),
                             ignore_index=True,
                             sort=False)

    # fna.to_csv("test.txt", sep="\t")
    results = pd.merge(results, fna)
    failed = results[(results.Success == 0) | (results.fna == False)]
    if(len(failed.index) > 0):
        print("Writing failed genomes file ({})".format(args.failed))
        failed.to_csv(args.failed, sep="\t", index=False)

    print("{} genomes downloaded.".format(str(sum(results.Success == 1))))
    print("{} genomes failed.".format(str(sum(results.Success == 0))))
    print("{} genomes have no .fna file.".format(str(sum(results.fna == False))))
