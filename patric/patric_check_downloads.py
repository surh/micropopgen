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
import os
import pandas as pd
import shutil


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Function that checks a folder that contains "
                          "genomes downloaded from PATRIC ftp site.")

    # Define required arguments
    required.add_argument("indir", help=("Directory containing patric "
                                         "downloads. Must contain a two-level "
                                         "directory structure where the first "
                                         "level is some arbitrary grouping "
                                         "and the second level contains one "
                                         "genome per directory."),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--outfile", help=("Name of file for output"),
                        type=str,
                        default="results.txt")
    parser.add_argument("--clean", help=("If included, remove failed "
                                         "folders"),
                        action="store_true",
                        default=False)

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    return args


def check_genomes_dirs(indir):
    """Takes a directory that contains a number of genome subdirectories,
    and checks that every genome subdirtectory has a .fna file

    Taken from patric_download_genomes.py"""

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

    specdirs = os.listdir(args.indir)

    Res = pd.DataFrame()
    for s in specdirs:
        s_dir = os.path.join(args.indir, s)
        res = check_genomes_dirs(s_dir)
        res['path'] = s_dir + '/' + res.ID
        Res = Res.append(res)

    # Write results
    Res.to_csv(args.outfile, sep="\t", index=False)

    if args.check:
        print("Cleaning")
        Res = Res[Res.fna != True]
        if Res.shape[0] > 0:
            for d in Res.path:
                print("Removing {}".format(d))
                shutil.rmtree(d)
        else:
            print("No directories to clean")
