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
import sqlite3
import argparse


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Take a table of eggNOG annotations and match them "
                          "with their description and COG categories. Table "
                          "must have a column called 'terms' with eggnog "
                          "annotations from emapper.py ")

    # Define required arguments
    required.add_argument("--input", help=("Input table file."),
                          required=True, type=str)
    required.add_argument("--db", help=("Path to sqlite3 db"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--outfile", help=("Outfiput file name."),
                        type=str,
                        default="eggnog_table.txt")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()

    # Read table
    Tab = pd.read_csv(args.input, sep="\t")

    # Get list of terms
    Tab['og'] = Tab.term.replace(to_replace="@\w*NOG",
                                 value="",
                                 regex=True)

    # Read db
    conn = sqlite3.connect(args.db)
    Dat = pd.read_sql_query("SELECT og, description, COG_categories FROM og",
                            conn)
    conn.close()

    # Subset data from db
    Dat = Dat.loc[Dat.og.isin(Tab.og), :]

    # Match table and DB
    Res = pd.merge(left=Tab, right=Dat,
                   how='left',
                   on='og',
                   sort=False).drop(columns='og')

    # Format COG cats
    Res.COG_categories.replace(to_replace="^\[u\'",
                               value="",
                               regex=True,
                               inplace=True)
    Res.COG_categories.replace(to_replace="'\]$",
                               value="",
                               regex=True,
                               inplace=True)

    # Write output file
    Res.to_csv(args.outfile, sep="\t")
