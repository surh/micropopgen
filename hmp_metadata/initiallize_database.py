#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

# The purpose of these script is to create 
# The skeleton of a database, and populate it with HMP data
import sqlite3
import argparse
import os

def create_metagenomes_database(path):
    
    # Create connection
    conn = sqlite3.connect(path)
    c = conn.cursor()
    
    # Create subject table
    create_table = '''CREATE TABLE subject
                    (subject_id INTEGER PRIMARY KEY,
                    gender TEXT,
                    internal_id TEXT,
                    disease_state TEXT,
                    population TEXT,
                    project TEXT)'''
    c.execute(create_table)
    
    # Create Sample table
    create_table = '''CREATE TABLE sample
                    (sample_id INTEGER PRIMARY KEY,
                    subject_id INTEGER,
                    body_site TEXT,
                    visit INTEGER,
                    age INTEGER,
                    SRS TEXT,
                    host_species TEXT,
                    replicate INTEGER,
                    FOREIGN KEY (subject_id)
                            REFERENCES subject(subject_id)
                    )
                    '''
    c.execute(create_table)
    
    # Create Run table
    create_table = '''CREATE TABLE run
                    (run_id TEXT PRIMARY KEY,
                    sample_id INTEGER,
                    subject_id INTEGER,
                    platform TEXT,
                    library_strategy TEXT,
                    SRR TEXT,
                    FOREIGN KEY (subject_id)
                        REFERENCES subject(subject_id),
                    FOREIGN KEY (sample_id)
                        REFERENCES sample(sample_id)
                    )'''
    c.execute(create_table)
    
    # Create tables
    conn.commit()
    conn.close()
    


if __name__ == '__main__':
    # Read arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--outpath", help = "Path of the file to create", type = str,
                          required = True)
    parser.add_argument("--overwrite", help = "If passed, it will overwrite any existing file",
                        action = "store_true")
    args = parser.parse_args()
    
    # Check that path exists
    head, tail = os.path.split(args.outpath)
    try:
        if not os.path.isdir(head):
            raise NotADirectoryError
    except:
        print("Path provided for output does not exist")
        raise
    
    # Check if file with same name exists
    try:
        fe = False
        if os.path.isfile(args.outpath):
            fe = True
        if fe and not args.overwrite:
            raise FileExistsError
        elif fe and args.overwrite:
            print("\tFile already existsts and it will be overwritten")
            os.unlink(args.outpath)
        elif not fe:
            print("\tFile does not exists and will be created")
    except:
        print("File already exists. Will do nothing")
        raise
    
    create_metagenomes_database(args.outpath)
    
