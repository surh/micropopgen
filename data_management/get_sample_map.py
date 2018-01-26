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


def read_samples(f):
    dat = pd.read_csv(f, header=None)
    samples = dat.loc[:, 0]
    return samples


def query_database(samples, db):
    # Create connection
    conn = sqlite3.connect(db)
    c = conn.cursor()

    # Create subject table
    create_table = '''SELECT (SRS, body_site) FROM sample
                    where SRS in ?'''
    res = c.execute(create_table, samples)

    conn.commit()
    conn.close()

    return res


if __name__ == "__main__":
    args = process_arguments()
    samples = read_samples(args.samples)
    print(samples)
    dat = query_database(samples, args.db)
    print(dat)
