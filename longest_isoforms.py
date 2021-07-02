#!/usr/bin/env python
# coding: utf-8
# created by Tdanis


import os
import sys
from itertools import groupby
import argparse
from collections import defaultdict

path_to_dir = '/home/thodoris/Desktop/prots/test'

fasta_filehandle = '../AED_less_1_mapped.fa' ######## this fasta file is a Maker annotation output , and especially after mapping to custom header proposed by Maker

purify_prots = path_to_dir+'/purify_prots.fasta'

def keep_uniq_prots(fasta_filehandle):          
   
    with open(fasta_filehandle, "r") as file, open('purify_prots.fasta', 'wt') as purify:
        
        dic_uniq_prots = defaultdict(list)
      
        for line in file:

            if line.startswith(">"):
                
                header = '\t'.join(line.split()[0:4:2])
                seq = next(file)
                dic_uniq_prots[seq.strip()].append(header[1:].split())
                
        purify_prot_ids = dict(map(lambda kv: (kv[0],sorted(kv[1], key = lambda item: float(item[1].split(':')[1]))), dic_uniq_prots.items())) ####  sort all values in list by AED score
        
        prots = 0
        for k,v in purify_prot_ids.items():

            aed = v[0][1].split(':')[1]
            
            if len(k) >= 30:  #### I choose to keep proteins above 30 AA 
                
                prots += 1
                
                purify.write('>'+v[0][0]+'\t'+str(len(k))+'\t'+str(aed)+'\n'+str(k)+'\n')
        print(f'Purified prots {prots}')
        return purify_prot_ids

uniq_prots = keep_uniq_prots(fasta_filehandle)



def keep_uniq_ids():
    
    # initiation of a dictinary for keeping longest isoforms
    
    dic = defaultdict(list) ### dict by gene id 

    c = 0
    for keys, values in uniq_prots.items():

        if len(values) >= 2:
#                 print(values)
            prot_ids = values[0][0]
            aed = values[0][1].split(':')[1]
            gene_ids = prot_ids.split('-')[0]

            dic[gene_ids].append((prot_ids,(len(keys),aed)))
        else:
#             print(values)
            aed = values[0][1].split(':')[1]
            prot_id = values[0][0]
            gene_ids = prot_id.split('-')[0]

            dic[gene_ids].append((prot_id,(len(keys),aed)))
    

    uniq_ids = []
    
##### iteration to keep the longest isoforms #######
    
    for keys, values in dic.items():


            values = sorted(values, key = lambda item: item[1][1], reverse = True)

            if len(values) < 2:
                
                prot = values[0][0]
                uniq_ids.append(prot)

            elif len(values) >= 2:

                prot_1 = values[0]
                prot_2 = values[1]

                if prot_1[1] == prot_2[1]:

                    prot1 = values[0][0]
                    uniq_ids.append(prot1)

                elif prot_1[1][1] > prot_2[1][1]:

                    prot2 = values[0][0]
                    uniq_ids.append(prot2)

                elif prot_1[1][0] < prot_2[1][0]:

                    prot3 = values[0][0]
                    uniq_ids.append(prot3)

                elif prot_1[1][1] == prot_2[1][1]:

                    prot4 = sorted(values, key = lambda item: item[1][0])[0][0]
                    uniq_ids.append(prot4)            
                else:
                    raise Exception ("Please check the type of the elemensr of these {values} ")

    return uniq_ids
    
uniq_ids = keep_uniq_ids()

##### retrieve back the prot sequencies and create a new fasta file
def main(purify_prots):
    
    with open(purify_prots, "r") as f, open('longest_isoforms.fasta', 'wt') as isoforms:
        
        gen = (x[1] for x in groupby(f, lambda line: line.startswith(">")))

        c = 0
        
        for header in gen:

            header = next(header)[1:].split()[0]
            seq = next(gen)

            if header in uniq_ids:
                c += 1

                isoforms.write(">"+header + '\n' + ''.join(list(seq)))
        print(f'Longest isoforms {c}')
main(purify_prots)
