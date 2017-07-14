#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

# Takes table indicating runs, and extracts all
# the metadata from the SRA API

import csv
import os
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
#sample_list_file = "/home/sur/micropopgen/data/HMP/test_hmiwgs.csv"
#sample_list_file = "/home/sur/micropopgen/data/HMP/HMIWGS_healthy.csv"
#sample_list_file = "/home/sur/micropopgen/data/HMP/test_hmp_cataloghe.csv"
sample_list_file = "/home/sur/micropopgen/data/HMP/hmp_catalogue_wgs_all_metadata.csv"
sample_col = 1
header = True
outdir = "out/"
total_runs_file = "total_runs.txt"
master_file = 'all_runs.txt'
skipped_file = "skipped_samples.txt"
keep_intermediate_files = True

# Prepare output
if not os.path.exists(outdir):
    os.mkdir(outdir)
else:
    print("Outdir ({}) already exists. Using it.".format(outdir))

sample_col -= 1
# print(sample_col)
with open(sample_list_file,'r') as infile:
    sample_colnames = None
    if header:
        sample_colnames = infile.readline()
    run_reader = csv.reader(infile, delimiter = ",")
    
    META = []
    RUNS = []
    SKIP = []
    for row in run_reader:
        sample = row[sample_col]
        print("Current sample is: {}".format(sample))
        
        # Prepare url for download and outfile 
        base_url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={}&result=read_run"
        search_url = base_url.format(sample)
        print(search_url)
        outfile = outdir + "/" + sample + ".csv"
        print(outfile)
        
        try:
            if os.path.exists(outfile):
                print("\tWARN: File {} exists already. Will not download".format(outfile))
            else:
                # Download and write result
                download = requests.get(search_url,timeout = 2)
                write_download(download, outfile)
                print("\tSucessfully downloaded...")
                
            # Get metadata from runs downloaded
            colnames, meta, nruns = process_ebi_metadata(outfile, sample)
            #print(len(meta))
            RUNS.append([sample,nruns])
            META.extend(meta)
            print("\tSaved metadata...")
        except (ConnectionError, Timeout):
            # Print warning if failed 
            print("\tWARN: Metadata for sample {} could not be downloaded ({}). Skipping.\n".format(sample,e))
            SKIP.append(row)
                
        # Clean
        if not keep_intermediate_files:
            os.unlink(outfile)

    # After iterating from all files, write output
    # Write total number of runs per sample
    write_table(outdir + "/" + total_runs_file, RUNS, header = ["Sample", "N.runs"], verbose = True)
    
    # Write master metadata
    write_table(outdir + "/" + master_file, META, header = colnames, verbose = True)
    
    # Write skipped files
    if len(SKIP) > 0:
        write_table(outdir + "/" + skipped_file, SKIP,  header = sample_colnames, verbose = True) 
