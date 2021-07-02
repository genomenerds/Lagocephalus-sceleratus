#!/usr/bin/env python
# coding: utf-8

# In[ ]:

# Created by: Teo Danis
# Descript: This script was used in order to keep the correrct orthogroups (for us) for the downstream analysis. It just needs the #Orthogroups.GeneCount.tsv file from Orthofinder and the Name of the column that you would like for process. For example our species is Lagocephalus #sceleratus and the column name is LSCEL.
#  

#===========================================================================================================



import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np

with open('Orthogroups.GeneCount.tsv') as f:

    df = pd.read_csv('Orthogroups.GeneCount.tsv', sep = '\t')


    new_df = df[(df[['LSCEL']] != 0).all(axis=1)] ###### keep rows !=0 for L. sceleratus column (choose your own species as it seems in headers)
    df1 = new_df[(new_df.iloc[:,1:-1] <= 1 ).all(axis=1)]   ####### keep all rows with no more than 1 gene 

    
    
    xi = []  
    
    total = 0
    
    species = []
    
    orthos = []
    
    
    
    list_of_species_genes = []
    
    for key, group in df1.groupby(['Total']):
        
        total += len(group)
        xi.append((key,len(group)))
        
        if key >= 26: ### keep orthologs that are represented at least from 26 per orthogroup
            
            orthos.append(len(group))

        l = list(group.iloc[:,1:-1].agg('sum'))
        
        species.append((np.count_nonzero(l),key,len(group)))
        names = list(group.columns[1:])
        
        tuples_species_genes = list(zip(names,l)) ## pairs of species and genes per orthogroup
        
        sorted_tuples_species_genes = sorted(tuples_species_genes, key = lambda tpl: float(tpl[1]), reverse = True)
        
        list_of_species_genes.append(sorted_tuples_species_genes)


    x, y = zip(*xi) # unpack a list of pairs into two tuples

    
    plt.rcParams["figure.figsize"] = (30,15)
    plt.bar(x,y, color = 'orange')
    
    ######################################################
    # add dominant species in every orthogroup
    ######################################################
    x_axis = 2
    for spe in list_of_species_genes:
        
        if 'LSCEL' in spe[1][0]:
            
            plt.text(x_axis-0.5, int(spe[0][1])+5, f'{spe[0][0], spe[0][1]}', rotation = 45, fontsize = 10)
            
        else:
            
            plt.text(x_axis-0.5, int(spe[0][1])+5, f'{spe[1][0], spe[1][1]}',  rotation = 45, fontsize = 10)
        
        x_axis += 1 
    ###################################################
    
    for sp in species:
        
        plt.text(sp[1]-0.5, sp[2]+1, f'{sp[0]}/34')
        
    plt.xlabel('Genes per orthogroup', fontsize= 20, labelpad = 15)
    plt.ylabel('Number of orthogroups', fontsize= 20, labelpad = 15)
    plt.show()
    print(f'* At the top of every bin there is one ratio which is the species that contribute to the orthogroups and the total number of them')
    print(f'Total number of orthogroups, {total}')
    print(f'Total number of orthogroups beyond 26 genes, {sum(orthos)}')
