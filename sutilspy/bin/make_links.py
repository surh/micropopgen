#!/usr/bin/env python
# Copyright (C) 2018 Sur Herrera Paredes

import argparse
import os
import sutilspy

def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Takes a file with a list of paths, "
                          "and creates a set of links to those paths "
                          "in the specified directory.")

    # Define required arguments
    required.add_argument("--paths", help=("A tab delimited file with "
                                           "the paths to be linked. "
                                           "it assumes first column has "
                                           "the paths, and no header."),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--col", help=("Column containing the paths to link."),
                        type=int,
                        default=1)
    parser.add_argument("--outdir", help=("Path where to write the links."),
                        type=str,
                        default='./')

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()

    # Get paths
    print("Reading list of paths")
    paths = sutilspy.io.return_column(args.paths, col=args.col, header=False)

    # Crete outdit
    print("Creating output dir")
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    print("Creating links")
    for f in paths:
        abspath = os.path.abspath(f)
        if not os.path.exists(abspath):
            raise FileNotFoundError("File {} not found".format(abspath))

        dst = '/'.join([args.outdir, os.path.basename(f)])

        os.symlink(src=abspath, dst=dst)
