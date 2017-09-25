
# coding: utf-8

# In[31]:


import pandas as pd
import argparse
import os
import sutilspy


# In[18]:


comparisons_file = '/home/sur/micropopgen/exp/2017/today4/comparisons.txt'
indir = '/godot/hmp/midas/merged.snps/'
map_file = '/home/sur/micropopgen/exp/2017/today4/full_map.txt'
mk_bin = '/home/sur/micropopgen/src/micropopgen/analyze_midas_snps/MKtest.py'
outdir = './out/'

# parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# required = parser.add_argument_group("Required arguments")
# required.add_argument("--comparisons", help = "File with comparisons to make",
#                       type = str, required = True)
# #     required.add_argument("--metadata_file", help = "Mapping file for samples", type = str,
# #                           required = True)
# #     required.add_argument("--group1", help = "Group1 of comparison", type = str,
# #                           required = True)
# #     required.add_argument("--group2", help = "Group2 of comparison",
# #                           required = True)
        
# #     parser.add_argument("--test", help = "Eventually specify test to perform",
# #                         default = "G", type = str)
# #     parser.add_argument("--outfile", help = "Output file with results",
# #                    default = "mk_results.txt", type = str)
# #     parser.add_argument("--min_count", help = "min depth at a position in a sample to consider that sample in that position",
# #                         default = 1, type = int)
# #     parser.add_argument("--nrows", help = "Number of gene positions to read",
# #                         default = float('inf'), type = float)
# #     parser.add_argument("--tables", help = "Output file for contingency tables",
# #                         default = "mk_tables.txt", type = str)
# #     parser.add_argument("--pseudocount", help = "Pseudocount value to use in contingency tables",
# #                         default = 1, type = int)
# args = parser.parse_args()
    


# In[19]:


comparisons = pd.read_csv(comparisons_file, sep="\t")
comparisons.head()


# In[30]:


for i, r in comparisons.iterrows():
    #print(i)
    #print(r)
    
    species_indir = ''.join([indir,'/',r['Species'], '/'])
    #print(species_dir)
    
    # Create output directory
    species_outdir = ''.join([outdir,'/',r['Species'],'/'])
    try:
        os.mkdir(species_outdir)
        print("Creating output directory")
    except FileExistsError:
        print("Directory already exists")
    
    suffix = r['A'] + '_' + r['B']
    suffix = suffix.replace(' ','.')
    #print(suffix)
        
    species_outfile = ''.join([species_outdir,
                               '/','mk_results.',
                               suffix,'.txt'])
    species_tables = ''.join([species_outdir,'/',
                              'mk_tables.',suffix,
                              '.txt'])
    species_log = ''.join([species_outdir,'/',
                           suffix,'.log'])
    species_err = ''.join([species_outdir,'/',
                           suffix,'.err'])
    
    cmd = [mk_bin, '--indir', species_dir,
          '--metadata_file', map_file,
           '--group1', r['A'],
           '--group2', r['B'],
           '--outfile', species_outfile,
           '--tables', species_tables,
           '1>', species_log,
           '2>', species_err]
    cmd = ' '.join(cmd)
    #print(cmd)
    sutilspy.io.run_command(cmd)
    
    #/micropopgen/src/micropopgen/analyze_midas_snps/MKtest.py --indir /godot/hmp/midas/merged.snps/Porphyromonas_sp_57899/ --metadata_file map.txt --group1 'Supragingival plaque' --group2 'Tongue dorsum' --outfile Porphyromonas_sp_57899_mk_results.txt --tables Porphyromonas_sp_57899_mk_tables.txt 1> log 2> err

    


# In[ ]:


i


