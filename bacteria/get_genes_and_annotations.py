# Copyright (C) 2017 Sur Herrera Paredes
import pandas as pd
import pybedtools as bed
from shutil import copyfile
import os
import argparse
# import numpy as np

def split_gene_annotations(functions, sep1=';',
                           sep2=':', sep3=',',
                           append_which=False):
    """Take annotation string from MIDAS database and expand it
    into a table of Annotation<->Gene"""

    # Split strings and create dictionary
    by_type = functions.split(sep=sep1)
    annotations = {}
    for t in by_type:
        db, annots = t.split(sep=sep2)
        annots = annots.split(sep=sep3)

        if append_which is True:
            annots = [sep2.join([db,i]) for i in annots]

        annotations[db] = annots

    # Get dataframe from dictionary
    d = _create_annotation_dataframe(annotations)
    return(d)

def _create_annotation_dataframe(annotations):
    """Expands annotation of one gene"""

    d = pd.DataFrame()
    for k,a in annotations.items():
        d2 = pd.DataFrame(data = a, columns=['Annotation'])
        d2['Type'] = k
        #print(d2)
        d = d.append(d2)
    return(d)


def get_fna(features, fasta, outdir, prefix):
    """Uses pybedtools to generate a fna file. A file with the nucleotide
    sequences of all CDS"""


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)

    # Define description
    parser.description = ("This script reads the MIDAS bacterial reference "
                          "database, extracts the gene annotations and "
                          "writes an fna file for all the genes\n\n"
                          "It must be noted that on version 1.2 of the "
                          "database the KEGG and EC numbers are simply "
                          "the GO numbers repeated, and those annotations "
                          "must not be used.\n\n"
                          "This script takes a directory as input, as well "
                          "as an instruction indicating if this directory "
                          "contains sub-directories for several genomes "
                          "or if the directory corresponds to a single "
                          "genome\n\n")

    # Define arguments
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--indir", help=("Path to directory where genome "
                                           "data is located. Each genome "
                                           "directory must contain a "
                                           "'features' file and an 'fna' "
                                           "file."),
                          required=True)
    required.add_argument("--type", help=("Indicates whether the directory "
                                          "is from a single genome or if "
                                          "it contains a number of genome "
                                          "sub-directories"),
                          choices=['single', 'multi'],
                          required=True)
    parser.add_argument("--which", help=("Indicates which annotation to "
                                         "extract. Eventually an 'all' "
                                         "option will be provided"),
                        choices=['GO', 'FIGFAM'],
                        default='GO')
    parser.add_argument('--outdir', help=("Directory where to save output"),
                        default="./", type=str)
    parser.add_argument("--overwrite", help=("If the output directory exists, "
                                             "passing this will overwrite any "
                                             "existing files with repeated "
                                             "names."),
                        action="store_true", default=False)
    parser.add_argument("--just_ID", help=("If passed the annotation ID "
                                           "allone will be returned. "
                                           "Alternativeley, the annotation "
                                           "database ID will be pre-pended "
                                           "to the annotation ID."),
                        default=False, action="store_true")
    parser.add_argument("--notdirs", help=("Indicates what to do in case "
                                           "that a 'multi' directory is "
                                           "passed, and it contains some "
                                           "non directory entries"),
                        choices=['ignore', 'fail'],
                        default='ignore')

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    # Prepare list of directories
    if args.type == 'multi':
        args.dirs = get_sample_dirs(args)
    elif args.type == 'single':
        args.dirs = [args.indir]
    else:
        print("ERROR: Incorrect type of directory passed")
        raise ValueError

    return args

def get_sample_dirs(args):
    """Gets list of subdirectories withina directory, and checks that
    there are no non-directory entries"""

    files = os.listdir(args.indir)

    dirs = []
    not_dirs = []
    for f in files:
        path = "".join([args.indir, "/", f])
        if os.path.isdir(path):
            dirs.append(path)
        else:
            not_dirs.append(path)

    if args.notdirs == 'fail':
        try:
            if len(not_dirs) > 0:
                raise ValueError
        except:
            print(("ERROR: The passed directory ({}) contains non-directory "
                   "entries").format(args.indir))
            raise

    return dirs



if __name__ == "__main__":
    """Main sscript body"""

    # Read arguments
    args = process_arguments()


    # Parameters
    infile = "/home/sur/micropopgen/exp/2017/today3/Zymomonas_mobilis_61858/genome.features"
    fasta = "/home/sur/micropopgen/exp/2017/today3/Zymomonas_mobilis_61858/genome.fna"
    which = "GO"
    outdir = "/home/sur/micropopgen/exp/2017/today3/out/"
    prefix = 'Zymomonas_mobilis_61858'
    append_which = True


    # In[ ]:


    # Read feature table
    os.mkdir(outdir)
    Feat = pd.read_csv(infile,sep="\t")
    Feat = Feat.head(n=10)
    # Feat.head()
    # split_gene_annotations(functions=Feat['functions'][0])


    # In[ ]:


    # Gett annotations
    ngenes = 0
    ncds = 0
    nannot = 0
    nwhich = 0
    Res = pd.DataFrame(columns=['Annotation','Type','Gene'])
    for i, r in Feat.iterrows():
        g = r['gene_id']
        a = r['functions']
        t = r['gene_type']
        ngenes = ngenes + 1
        # print(a)
        # print("============")

        # Keep only CDS
        if t != 'CDS':
            continue
        ncds = ncds + 1
        # Skip unannotated genes
        if pd.isnull(a):
            continue
        nannot = nannot + 1

        d = split_gene_annotations(functions=a,
                                   append_which=append_which)
        d['Gene'] = g

        # Select annotation
        if not pd.isnull(which):
            d = d.loc[ d.Type == which, :]

        # Append only if it has rows
        if len(d.index) > 0:
            nwhich = nwhich + 1
            Res = Res.append(d)

    Res = Res.drop(['Type'], axis = 1)



    # In[ ]:


    print(ngenes,ncds,nannot,nwhich)


    # In[ ]:


    Res


    # In[ ]:


    # Get BED format dataframe
    Bed = Feat[ ['scaffold_id', 'start', 'end', 'gene_id', 'gene_type', 'strand']].copy()
    Bed.reset_index(drop=True)
    print(Bed)
    Bed.loc[ Bed.strand == '+', 'start'] = Bed.loc[ Bed.strand == '+', 'start'] - 1
    Bed.loc[ Bed.strand == '-', 'start'] = Bed.loc[ Bed.strand == '-', 'start'] - 1
    Bed = Bed.loc[ Bed.gene_type == 'CDS', ]
    Bed


    # In[ ]:


    (Bed.end - Bed.start)/3


    # In[ ]:


    # create BedTool and obtain sequences
    Bed = bed.BedTool.from_dataframe(Bed)
    Bed.sequence(fi=fasta, s=True, name=True)


    # In[ ]:


    # Show results
    print(open(Bed.seqfn).read())


    # In[ ]:


    Bed.seqfn


    # In[ ]:


    outfile = ''.join([outdir,'/',prefix,'.cds.fna'])


    # In[ ]:


    copyfile(src=Bed.seqfn,dst=outfile)
