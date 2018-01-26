#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

import sqlite3
import argparse
import pandas as pd


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Takes list of samples and database")

    # Define required arguments
    required.add_argument("--samples", help=("File containing sample IDs, "
                                             "one per line"),
                          required=True, type=str)
    required.add_argument("--db", help=("sqlite3 database"))

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


def read_samples(file):
    samples = []
    with open(file) as f:
        for l in f.readline():
            samples.append(l)

    return samples


def query_database(args)

if __name__ == "__main__":
    args = process_arguments()
    samples = read_samples(args.samples)
    dat = query_database(samples, args.db)
