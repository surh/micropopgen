
# coding: utf-8

# In[ ]:


# Copyright (C) 2017 Sur Herrera Paredes
import pandas as pd
import pybedtools as bed
from shutil import copyfile
import os
# import numpy as np


# In[ ]:


def split_gene_annotations(functions, sep1=';',
                           sep2=':', sep3=',',
                           append_which = False):
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


# In[ ]:


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

