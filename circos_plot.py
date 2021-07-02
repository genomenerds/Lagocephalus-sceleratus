
                                                   ##############Circos Plot############


################# required files 

#file_1 = LSCEL__v__TETNG.tsv  -----> orthofinder results
#file_2 = nigro_coordinates.txt  -----> biomart files
#file_3 = map  -----> Maker map file (follow basic tutorial)
#file_4 = TETNG.fa -----> uniprot proteome
#file_5 = lscel.gff3 -----> annotation file (Maker2 gff3)


#resulting_file1 = first_part_pands.txt
#resulting_file2 = finals_bands.txt




######### get one-to-one orthologs from orthofinder tsv results ############

def count_one2one_orthologs():
    
    with open("LSCEL__v__TETNG.tsv") as tsv1:
        os1 = 0
        putative_genes = []

        for line1 in tsv1:
            line1 = line1.strip()
            line1 = line1.replace(',', '')
            line1 = line1.split()
            if line1[0].startswith('Orth'):
                h1 = line1[2]

            if len(line1)== 3 and line1[0].startswith('OG'):
                os1 += 1
                putative_genes.append((line1[1],line1[2]))

        return putative_genes
    
final_orthologs = dict(count_one2one_orthologs()) ##### resulting dictionary, for example  ( Lscel_00013489-RA': 'TETNG16956') 


################## creating dictionary with the coordinates for the potential subject species ######### Chromosome : (start, end)


with open('nigro_coordinates.txt') as fasta:  ######## chromosome coordinates for T. nigroviridis retrieved from Biomart , Ensembl version 98 (Gene_Ensembl_id, chromosome, start, end)
    
    nigro_dict_coords = {line.strip().split()[0]:(line.strip().split()[1],line.strip().split()[2],line.strip().split()[3]) for line in fasta if line.strip().split()[1].isdigit() or  line.strip().split()[1] == 'Un_random'}




######## 41 largest contigs which represent 91% of our genome ########

    contigs = ['scaffold_101:1.0-17057981.0_pilon', 'scaffold_106:1.0-15598540.0_pilon', 'contig_127:1.0-14919945.0_pilon', 'contig_22:4.0-14575672.0_pilon', \
               'contig_16:1.0-13974425.0_pilon', 'contig_118:1.0-13832250.0_pilon', 'contig_174:1.0-13151330.0_pilon', 'contig_533:1.0-12440543.0_pilon', 'scaffold_67:1.0-12279225.0_pilon', \
               'contig_8:1.0-12122666.0_pilon', 'contig_91:1.0-11773260.0_pilon', 'scaffold_38:1.0-11588445.0_pilon', 'contig_90:1.0-11438290.0_pilon', 'contig_93:1.0-11311637.0_pilon', \
               'contig_123:1.0-11274211.0_pilon', 'contig_9:1.0-11125835.0_pilon', 'contig_15:1.0-11122009.0_pilon', 'contig_427:384.0-9460825.0_pilon', 'contig_19:1.0-9019598.0_pilon', \
               'contig_599:6529.0-8758774.0_pilon', 'contig_143:1.0-8481133.0_pilon', 'contig_2:1.0-8206783.0_pilon', 'contig_28:1.0-8182740.0_pilon', 'scaffold_125:1.0-7722920.0_pilon', \
               'scaffold_838:1.0-6380105.0_pilon', 'contig_122:1.0-6309019.0_pilon', 'contig_94:1.0-4884195.0_pilon', 'contig_14:1.0-4861701.0_pilon', 'contig_124:1.0-4492329.0_pilon', \
               'scaffold_115:1.0-3763395.0_pilon', 'contig_20:1.0-3597822.0_pilon', 'contig_6:1.0-3057726.0_pilon', 'contig_3:1.0-3047716.0_pilon', 'contig_4:1.0-2449395.0_pilon', \
               'contig_153:2264884.0-4640337.0_pilon', 'contig_153:1.0-2262385.0_pilon', 'contig_200:1.0-2164460.0_pilon', 'contig_23:1.0-2122701.0_pilon', 'contig_191:1.0-1908425.0_pilon', \
               'scaffold_1001:1.0-1813021.0_pilon', 'contig_142:1.0-1794019.0_pilon']


##################################################33




################### getting the gene for the 41 largest contigs from L. sceleratus, using map file from maker


def get_genes_for_largest_contigs():
    with open("/home/thodoris/Desktop/orthologs_40_contigs/map") as map_file:

        c = 0
        contig_code = []
        for line in map_file:
            c+= 1 
            line = line.strip()
            contig = (line.split()[0]).split('-')[1:3]
            code = line.split()[1]
            tpl = tuple(('-'.join(contig), code))
            if tpl[0] in contigs:
                contig_code.append((tpl[0], code))
        return contig_code
    
contigs_codes_41_contigs = get_genes_for_largest_contigs()



res = dict((k,v) for k,v in final_orthologs.items()) ####### dictionary with L. sceleratus and T. nigroviridis genes ({'Lscel_00013489-RA': 'TETNG16956'.....})




# finalGenes = {}  ## nigro, lscel
# for contig, code in contigs_codes_41_contigs:
# #     print(contig,code)
#     if code in res.keys():
#         finalGenes[res[code]] = code


finalGenes = {res[code]:code for contig, code in contigs_codes_41_contigs if code in res.keys()}  ############## inversion of res dictionary 
        
print(len(finalGenes)) 

print(list(finalGenes.items())[:10])
        
ids_ensembl_nigro = {}

with open('TETNG.fa') as fasta:
    
    for line in fasta:

        line = line.strip()

        if line.startswith('>'):

            line = line.split('|')
            key = line[0].replace('>', '').replace(" ", "")
            value = line[1].replace(" ", "")
            ids_ensembl_nigro[key] = value

    
    
########################    pairing L. sceleratus genes with T. nigroviridis genes

final_pairs_lago_ensIDS = {lago:ids_ensembl_nigro[nigro] for nigro, lago in finalGenes.items() if nigro in ids_ensembl_nigro.keys()}

#################
# print(final_pairs_lago_ensIDS)





############## getting coordinates for L. sceleratus genes 


def get_coordinates_orthos():
    
    with open('Lscel.gff3') as gff:
        c = 0
        coordinates = {}

        for line in gff:

            if line.startswith('#'): continue
            line = line.strip().split('\t')

            if 'mRNA' in line:

                iD = line[-1].split(';')[0][3:]
                coordinates[iD] = (line[3], line[4])

        return coordinates

lago_dict_coords = get_coordinates_orthos()


lago_code_contig={v:k for k,v in contigs_codes_41_contigs}



############# constructing the first part of circos input, cords and chromosomes

with open('first_part_bands.txt', 'wt') as bands:
    
   
    for k,v in final_pairs_lago_ensIDS.items():
        
        if k in finalGenes.values():
            lago_id = k
            lago = lago_code_contig[lago_id].split('.')[0].replace(":","_")
            lago_start = lago_dict_coords[lago_id][0]
            lago_end = lago_dict_coords[lago_id][1]
            
            bands.write(f'{lago} {lago_start} {lago_end} {v}\n')
            #print(f'{lago} {lago_start} {lago_end} {v}\n')


print("first part of bands done!")  
####################################


################# finalizing circos bands, checking of the pairs contain of not uplaced regions which must be discarded


with open('first_part_bands.txt') as first_bands, open('final_bands.txt', 'wt') as final_bands:
   
    discard = []
    for line in first_bands:
        
        line = line.strip().split()
        
        nigro_id = line[-1]
        
        if  nigro_id in nigro_dict_coords.keys():
                  
            lago = line[0]
            lago_start = line[1]
            lago_end = line[2]
            
            nigro = '_'.join(('nigro',nigro_dict_coords[nigro_id][0]))
            nigro_start = nigro_dict_coords[nigro_id][1]
            nigro_end = nigro_dict_coords[nigro_id][2]
            
            final_bands.write(f'{lago} {lago_start} {lago_end} {nigro} {nigro_start} {nigro_end}\n')

        else :
            
            discard.append(lago_id) #### all the unplaced chrs

print("final part of bands done!")





# print(len(final_pairs_ensIDS_lago))
# print(list(final_pairs_ensIDS_lago.items())[:10])
# print(list(lago_dict_coords.items())[:10])
# print(len(finalGenes))
# print(list(finalGenes.items())[:10])
# print(nigro_dict_coords)
