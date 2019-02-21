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

import pandas as pd
import argparse
import math


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Takes a single file containing a table, and "
                          "splits it into multiple files.")

    # Define required arguments
    required.add_argument("FILE", help=("Tab-delimited file, and the first "
                                        "row must contain the headers"),
                          type=str, nargs=1)

    # Optional arguments
    parser.add_argument("--nchunks", help=("Number of new files to create"),
                        type=int, default=0)
    parser.add_argument("--nlines", help=("Number of lines of each new file"),
                        type=int, default=0)
    parser.add_argument("--outdir", help="Output directory",
                        type=str, default="./")
    parser.add_argument("--prefix", help=("Filename to write the "
                                          "concatenated results"),
                        type=str,
                        default="tab")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    if args.nchunks > 0 and args.nlines > 0:
        raise ValueError("ERROR: only one of --nchunks and --nlines must be "
                         "defined.")
    if args.nchunks == 0 and args.nlines == 0:
        raise ValueError("ERROR: one of --nchunks and --nlines must be "
                         "defined")

    return args


def split_table(Tab, step_size, outdir, prefix):
    nrows = Tab.shape[0]
    chunk = 1
    for i in range(0, nrows, step_size):
        tab = Tab.iloc[i:(i+step_size), :]
        filename = ''.join([outdir, "/",
                            prefix,
                            str(chunk),
                            '.txt'])
        print("Writing chunk {}".format(str(chunk)))
        tab.to_csv(filename, sep="\t", index=False)
        chunk = chunk + 1

    return


if __name__ == "__main__":
    args = process_arguments()

    Tab = pd.read_csv(args.FILE, sep="\t")
    nrows = Tab.shape[0]
    if args.nchunks > 0:
        if args.nchunks > nrows:
            args.nchunks = nrows
        step_size = nrows / args.nchunks
        step_site = math.ceil(step_size)

        split_table(Tab, step_site, args.outdir, args.prefix)

    elif args.nlines > 0:
        split_table(Tab, args.nlines, args.outdir, args.prefix)
