#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

# Takes table indicating runs, and extracts all
# the metadata from the SRA API

import csv
import wget
import os
# from wget import bar_adaptive
import requests
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
                
                res.append(row)
        meta_file.close()
    return(header,res,nruns)
                    
def write_table(outfile,rows, header = None, delimiter = "\t", verbose = False):
    with open(outfile,'w') as out_fh:
        writer = csv.writer(out_fh,delimiter = '\t')
        if verbose:
            print("\tWriting {}".format(outfile))
            
        nlines = 0
        if header is not None:
            writer.writerow(header)
            nlines += 1
        for row in rows:
            writer.writerow(row)
            nlines += 1
    out_fh.close()

    if verbose:
        print("\t\tWrote {} lines".format(nlines))
    
    return(nlines)        
                
def write_download(download,outfile):
    with open(outfile,'w') as out_fh:
        out_fh.write(download.text)
    out_fh.close()
    
    

###################################


# Global variables
# Eventually command line parameters
run_list_file = "/home/sur/micropopgen/data/HMP/test_hmiwgs.csv"
#run_list_file = "/home/sur/micropopgen/data/HMP/HMIWGS_healthy.csv"
sample_col = 1
header = True
outdir = "out/"
total_runs_file = "total_runs.txt"
master_file = 'all_runs.txt'
keep_intermediate_files = True

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
    
    META = []
    RUNS = []
    SKIPPED = []
    for row in run_reader:
        sample = row[sample_col]
        print("Current sample is: {}".format(sample))
        
        # Download run list and metadata for sample 
        base_url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={}&result=read_run"
        search_url = base_url.format(sample)
        
        #base_url = 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term='
        #search_url = base_url + sample
        
        print(search_url)
        outfile = outdir + "/" + sample + ".csv"
        if os.path.exists(outfile):
            print(outfile)
            print("\tWARN: File {} exists already. Will not download".format(outfile))
        else:
            #meta_file = wget.download(search_url, outfile, bar_adaptive)
            try:
                download = requests.get(search_url,timeout = 2)
                write_download(download, outfile)
                
                # Get metadata from runs downloaded
                header, meta, nruns = process_ebi_metadata(outfile, sample)
                #print(len(meta))
                RUNS.append([sample,nruns])
                META.extend(meta)
            except (ConnectionError, Timeout):
                print("WARN: Metadata for sample {} could not be downloaded ({}). Skipping.\n".format(sample,e))
                SKIP.append(sample)
            
        #print(meta_file)
        
        # Clean
        if not keep_intermediate_files:
            os.unlink(outfile)

    #print(len(header))
    #print(len(META))
    # Write total number of runs per sample
    write_table(outdir + "/" + total_runs_file, RUNS, header = ["Sample", "N.runs"], verbose = True)
    
    # Write master metadata
    write_table(outdir + "/" + master_file, META, header = header, verbose = True)    

