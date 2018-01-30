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
    required.add_argument("--outfile", help=("Outfile name"),
                          type=str, default='map.txt')

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


def read_samples(f):
    dat = pd.read_csv(f, header=None)
    samples = dat.loc[:, 0]
    samples = [tuple(x) for x in [samples]]
    return samples[0]


def query_database(samples, db):
    # Create connection
    conn = sqlite3.connect(db)
    c = conn.cursor()

    # Create subject table
    dat = []
    notfound = []
    for s in samples:
        query_db = 'SELECT SRS, body_site FROM sample WHERE SRS=?'
        c.execute(query_db, (s,))
        res = c.fetchall()
        if(len(res) > 0):
            dat.append(res)
        else:
            notfound.append([s])

        res = []

    return dat, notfound


def write_output(dat, outfile):
    with open(outfile, 'w') as o:
        for l in dat:
            # print(l)
            s = "\t".join(l[0]) + "\n"
            # print(s)
            o.write(s)


if __name__ == "__main__":
    args = process_arguments()
    samples = read_samples(args.samples)
    # print(samples)
    dat, notfound = query_database(samples, args.db)
    # print(dat)
    print(notfound)
    write_output(dat, args.outfile)
