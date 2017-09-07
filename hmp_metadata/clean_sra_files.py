#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

import os
import sutilspy

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    # Required arguments
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--indir","-i", help = "Directory containing .sra files from all the runs to be processed",
                          type = str, required = True)
    required.add_argument("--map","-m", help = "Input tab delimited file that has list of SRA IDs to be eliminated",
                        type = str, required = True)

    # Optional arguments
    parser.add_argument("--col", help = "Column where the SRA ID is located in map", type = int,
                        default = 1)
    parser.add_argument("--header", help = "Flag indicating whether table has headers in the first row",
                        action = "store_true")
    parser.add_argument("--failed", help = "File to store the failed samples", type = str, default = 'failed.txt')

    args = parser.parse_args()

    sra_ids = sutilspy.io.return_column(args.map,args.col,header = args.header)

    for sra in sra_ids:
        sra_file = args.indir + '/' + sra + '.sra'
        print("Removing {}".format(sra_file))
        os.remove(sra_file)
