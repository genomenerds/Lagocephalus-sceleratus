#!/usr/bin/env python
# coding: utf-8

# In[ ]:

# create by Tdanis

from collections import defaultdict, Counter
import itertools
import copy
import pandas as pd




with open('dump.blast_output.mci.I30') as file, open('unfiltered.txt', 'wt') as txt:
    
    
    species_ids = "0 1 2 3 4 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22 24 25 26 27 28 29 31 32 33 34" ## provide species ids
    
    sps_header = "\t".join(species_ids.split())

    sp_dict = {sp:0 for sp in species_ids.split()}

    mega_list = []
    mega_dict ={}
    c = 0
    
    len_lines = []
    for line in file:

        c += 1
        line = line.strip()
        len_lines.append(len(line.split()))
        for sp_id in line.split():

            if sp_id.split('_')[0] in sp_dict.keys():
                          
                spID = sp_id.split('_')[0]        
                
                sp_dict[spID] += 1
                sp_dict.copy()
#         print(sp_dict)
        species_dict = copy.deepcopy(sp_dict)
        mega_list.append(species_dict)
        sp_dict = {sp:0 for sp in species_ids.split()}

    
    header = 'Desc'+'\t'+'Family ID'+'\t'+'\t'.join(species_ids.split())
    txt.write(header+'\n')
    
    for l,p in zip(mega_list, len_lines):
    
        assert sum(l.values()) == p , "Oh no! This assertion failed!. The number of genes are not the equivalent with the those of unfilterder file"
    
    
    family_id = 0
    for dict_with_genes in mega_list:
        family_id += 1
        line = '\t'.join(str(gene) for gene in dict_with_genes.values())
        txt.write('(null)' + '\t' + str(family_id) + '\t'+ line+'\n')
        
print('Done filtering')


with open('species.txt') as sp:   #### provide a species txt tab separated 
    
    sps = [x for line in sp for x in line.strip().split()]

df = pd.read_csv('unfiltered.txt', sep='\t')
df
df1 = df[(df.iloc[:,2:] <= 100 ).all(axis=1)] ### below 100 genes per family for every speceis
df1.columns = sps

df2 = df[~(df.iloc[:,2:] <= 100 ).all(axis=1)] ### over 100 genes for at leat one family

with open('filtered.txt', 'wt') as txt_f:
    txt_f.write(df1.to_csv(sep='\t', index = False, header=True))
    
with open('filtered_over100.txt', 'wt') as txt_f:
    txt_f.write(df2.to_csv(sep='\t', index = False, header=True))
    
print('Done families seperation')
