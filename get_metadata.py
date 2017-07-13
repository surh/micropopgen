#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

# Takes table indicating runs, and extracts all
# the metadata from the SRA API

import csv
import wget
import os
from wget import bar_adaptive
import re

################ Modules ##########
def process_ebi_metadata(infile,accession,accession_col = 4):
    accession_col -= 1
    with open(infile,'r') as meta_file:
        reader = csv.reader(meta_file,delimiter = '\t')
        header = next(reader)
        header.append('subject_id')
        print("\tSeaching for accession {} in column {}".format(accession,header[accession_col]))
        #print(header[accession_col])
        nruns = 0
        res = []
        for row in reader:
            if row[accession_col] == accession:
                nruns += 1
                #print(row[22]) 
                
                # Search for subject id
                title = row[22]
                m = re.search('containing sample (\d+) from participant (\d+)', title)
                if m is not None:
                    #print(m.group(1))
                    row.append(m.group(1))
                else:
                    row.append('NA')
        meta_file.close()
    return([header,res,nruns])
                    


###################################


# Global variables
# Eventually command line parameters
run_list_file = "/home/sur/micropopgen/data/HMP/test_hmiwgs.csv"
sample_col = 1
header = True
outdir = "out/"

# Prepare output
if not os.path.exists(outdir):
    os.mkdir(outdir)
else:
    print("Outdir ({}) already exists. Using it.".format(outdir))

sample_col -= 1
# print(sample_col)
with open(run_list_file,'r') as infile:
    if header:
        infile.readline()
    run_reader = csv.reader(infile, delimiter = ",")
    for row in run_reader:
        sample = row[sample_col]
        print("Current sample is: {}".format(sample))
        
        base_url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={}&result=read_run"
        search_url = base_url.format(sample)
        
        #base_url = 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term='
        #search_url = base_url + sample
        
        print(search_url)
        outfile = outdir + "/" + sample + ".csv"
        print(outfile)
        #meta_file = wget.download(search_url, outfile, bar_adaptive)
        #print(meta_file)
        process_ebi_metadata(outfile, sample)
        



