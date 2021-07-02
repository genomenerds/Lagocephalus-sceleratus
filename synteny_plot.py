import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib import collections  as mc

# The commands used for last and get the tab file
# Example of alignment for L. sceleratus and T. nigroviridis (as reference)

#lastdb -P0 -uMAM4 -R01 -v Lago-MAM4 draft.fa && \
#last-train -P0 --revsym --matsym --gapsym -E0.05 -C2 Lago-MAM4 Tnigro.fasta > Lago-MAM4-nigro.mat && \
#lastal -m100 -P0 -E0.05 -C2 -p Lago-MAM4-nigro.mat Lago-MAM4 Tnigro.fasta | last-split -m1 > Lago-nigro-1.maf && \
#maf-swap Lago-nigro-1.maf | awk '/^s/ {$2 = (++s % 2 ? "nigro." : "Lago.") $2} 1' | last-split -m1 | maf-swap > Lago-nigro-2.maf && \
#last-postmask Lago-nigro-2.maf | /mnt/big/WholeGenomeAlignment/last-1047/scripts/maf-convert -n tab | awk -F'=' '$2 <= 1e-3' > Lago-nigro.tab 

tab_file = 'path/to/tab/file' # output from last 
txt_file = 'path/to/txt/file' # headers from fasta file of reference genome (grep it)
reference_genome = 'path/to/reference/genome/reference.fasta'
draft_assembly = 'path/to/your/draft/assembly/draft.fasta'
tsv_file = 'path/to/output/tsv/file' 
assembly_blocks_txt = 'path/to/blocks/output/file'




def parsing_txt(txt_file):
   
    # Dictionary with ids as keys and corresponding chromosome as value for reference genome
    with open(txt_file, "r+") as file:
        pairs_list = []
        for line in file:
            
            if "chromosome" in line.split(" "):
                
                line = line.split(" ")
                indx = line.index("chromosome")
                
                indx2 = indx + 1 
                Chr = line[indx]
               
                num = line[indx2].split(",")
   
                final_chrom = Chr + "_" + num[-2]
                pairs_list.append((line[0][1:], final_chrom))

        dict1=dict(pairs_list)
          
        return dict1
           



def parsing_maftab(cutoff,tab_file,assembly_blocks_txt):
    
    # Creation of the blocks file
    with open(tab_file, "r+") as f, open(assembly_blocks_txt, "wt") as blocks_file:
        
        blocks_file.write(f'chrom\tcontig\tchrom_len\tcontig_len\tchrom_start\tchrom_end\tcontig_start\tcontig_end\tmismap\n')

        for line in f:
           
            if line.startswith('#'): continue
            line  =line.strip().split()
            contig = line[1][9:] ### search Lago.contig_1006:1.0-8269.0_pilon without the prefix Lago.
            contig_start = line[2]
            contig_alignmentSpan = line[2]
            contig_length = line[5]
             
            chrom = line[6][5:] ## without the prefix fugu

            if chrom in ids_pairs.keys():
                
                chrom = ids_pairs[chrom]
                chr_start = line[7]
                chr_alignmentSpan = line[8]
                chrom_length = line[10]
                mismap = float(line[-1].split('=')[-1])
                
                if mismap <= cutoff:

                    blocks_file.write(f'{chrom}\t{contig}\t{chrom_length}\t{contig_length}\t{chr_start}\t{int(chr_start)+int(chr_alignmentSpan)}\t{contig_start}\t{int(contig_start)+int(contig_alignmentSpan)}\t{mismap}\n')




def count_bases():
    
    # Returns:
    # 1. Dict with total bases per contig/scaffold of draft assembly
    # 2. Total length of assembly
    with open(draft_assembly, "r+") as f2:
        
        dic = {}
        
        for line in f2:
            line = line.strip()
            if line.startswith(">"):
                line = line.split(" ")
                header = line[0][1:] 
                dic[header] = 0
            
            else:
                dic[header] += len(line)
                 
        total = sum(dic.values())
        sorted_reversed = dict(sorted(dic.items(), key=lambda item: item[1], reverse =True))
        total = sum(sorted_reversed.values())
        
        return sorted_reversed, total




def count_bases_Reference(tsv_file):
    
    # Returns:
    # 1. Dict with total bases per contig/scaffold of draft assembly
    # 2. Total length of assembly
    # 3. Count chromosomes
    # Also, creates a tsv file with total length of each chromosome in reference genome (in Mb)
    with open(reference_genome, "r+") as f2, open(tsv_file, 'wt') as tsv:
        
        dic2 = {}
        for line in f2:
            
            line = line.strip()
            
            
            if line.startswith(">") and "chromosome" in line:

                line = line.split(" ")

                header = "_".join((line[3:5])).replace(",","")

                dic2[header] = 0

            elif line.startswith(">") and "chromosome" not in line:
                break

            else:
                dic2[header] += len(line) / 10**6


        total = sum(dic2.values())
        
        num_chrs = {chrom.split('_')[-1]:'chromosome_'+chrom.split('_')[-1] for chrom, length in dic2.items() if len(chrom.split('_')) <= 2 and chrom.split('_')[-1].isdigit()}
        
        tsv.write("chr" + "\t" + "size(Mb)" + '\n')
        for ch, length in dic2.items():

            tsv.write(ch + '\t' + f"{length:.2f}" + '\n')
     
        return dic2 , total, num_chrs



# Call the above functions
ids_pairs = parsing_txt(txt_file)
parsing_maftab(1e-03,tab_file,assembly_blocks_txt)
dict_contig_len, total_bases = count_bases()
bases_for_chrs, total, d3 = count_bases_Reference(tsv_file)




# Plot construction 
figsize = (12, 14)
fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)
    
longest_41_draft_contigs = dict(list(dict_contig_len.items())[:41]) ## Plot the longest 41 contigs from draft to reference
print(longest_41_draft_contigs)


chrs_contigs_regions = pd.read_csv(assembly_blocks_txt, sep='\t', header=0)
chrs_contigs_regions = chrs_contigs_regions.drop(["chrom_len", 'mismap', 'contig_len', 'contig_end', 'contig_start', 'contig_len' ], axis = 1)
chrs_contigs_regions['width'] = chrs_contigs_regions.chrom_end - chrs_contigs_regions.chrom_start


chrs_contigs_regions  = chrs_contigs_regions[chrs_contigs_regions.contig.apply(lambda x: x in longest_41_draft_contigs.keys() )]



colors = ['gray', 'silver', 'rosybrown', 'firebrick','red','maroon', \
          'chocolate', 'darksalmon', 'sienna','indianred','sandybrown', \
          'tan', 'yellow', 'gold','peachpuff', 'darkorange','darkgoldenrod',\
          'khaki','darkkhaki', 'chartreuse','olivedrab', 'mediumspringgreen',\
          'darkgreen', 'darkslategray', 'darkcyan','c', 'aqua', 'deepskyblue',\
          'lightsteelblue', 'cornflowerblue', 'lavender', 'navy','blue', \
          'mediumpurple','darkorchid','plum', 'mediumvioletred', 'palevioletred', \
          'purple', 'fuchsia', 'indigo' ]
                               
contigs_list = list(longest_41_draft_contigs.keys())

contigs_colors = dict(zip(contigs_list,colors))

chrs_contigs_regions['colors'] = chrs_contigs_regions['contig'].apply(lambda x: contigs_colors[x])


    

refChrs_lengths = pd.read_csv(tsv_file, sep='\t', header=0)  
refChrs_lengths['start'] = int(1)

start_end_limits = list(zip(refChrs_lengths['start'], refChrs_lengths['size(Mb)']))


sorted_list = sorted(map(lambda x : x.split("_"),list(set(chrs_contigs_regions['chrom']))),  key = lambda x: int(x[1]))


yvalue = 0
chrom_ypos = {}
chrom_centers = {}
start_yposition = {}



for chrom in sorted_list:
    chrom="_".join((chrom[0],chrom[1]))
    chrom_centers[chrom] = yvalue + 0.5
    start_yposition[chrom] = yvalue + 2 
    yvalue += 2
    



def contiguity_plot():

    for chrom, group in chrs_contigs_regions.groupby('chrom'):
        
        yrange = (chrom_centers[chrom], 0.5)
        group['width'] = group['width'].apply(lambda x: x)
        xranges = group[['chrom_start', 'width']].values

        yield BrokenBarHCollection(xranges, yrange, facecolors=group['colors'])



for collection in contiguity_plot():

    ax.add_collection(collection)

#Set appropriate ticks for x axis
ls = list(refChrs_lengths['size(Mb)'].apply(lambda x: x*(10**6)))
prevList = ls[0]

new_points = {}

a = [((ls[i])-ls[i+1]) for i in range(0, len(ls) -1 )]

for i in range(0, len(ls) -1 ):
    
    diff = abs(ls[i] - ls[i+1])
    if diff > 2*10**6:
        
        new_points[ls[i]] = 2 
        
    else:
        
        new_points[ls[i]] = 15


# Set yticks and labels
points = refChrs_lengths['size(Mb)'].apply(lambda x: x*(10**6))
labels = [str(point/10**6) for point in new_points.keys()]


points_labels = dict(zip(points,labels))

lengths = [new_points[k] for k,v in points_labels.items()]

chr_end_limits = dict(list(zip(refChrs_lengths['chr'], refChrs_lengths['size(Mb)']*(10**6))))

ax.set_yticks([chrom_centers[i] for i in list(set(chrs_contigs_regions['chrom']))])
ax.set_yticklabels(list(set(chrs_contigs_regions['chrom'])))



# Set lines at the end of every chromosome
y0 = 0.5
ytext = 1
y1 = 1
lines = []
li = []

for ch in sorted_list:

    ch="_".join((ch[0],ch[1]))
    
    if ch in chr_end_limits.keys():
        
    
        x0 = chr_end_limits.get(ch)
             
        x1 = chr_end_limits.get(ch)
        
        plt.text(x0, ytext, 'end')
        
        lines.append([(x0,y0),(x1,y1)])
        li.append([(ch,(x0,y0),(x1,y1))])
        
        y0 += 2
        ytext += 2
        y1 += 2
                   

lc = mc.LineCollection(lines, colors='black', linewidths=2)
ax.add_collection(lc)

# Axis 
ax.axis('tight')
ax.set_ylabel('Reference chromosomes', labelpad = 18, fontsize = 22)
ax.set_xlabel('Bases in Gbases',labelpad = 18, fontsize = 22)
ax.set_title('Synteny plot', {'fontsize': 20, 'fontweight' : 'bold' }, pad = 18)
markers = [plt.Line2D([0,0],[0,0],color=color, marker='o', linestyle='') for color in contigs_colors.values()]
plt.legend(markers, contigs_colors.keys(), numpoints=1, bbox_to_anchor=(1.6,1.0), loc='upper right')
plt.show()

# If you want to get the results of the plot in a file
#chrs_contigs_regions.to_csv('path/to/csv/file/results.tsv', sep='\t')
