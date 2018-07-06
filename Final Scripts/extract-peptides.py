import os
import csv
from collections import Counter

organisms = {}

dstDir = '../data/Uniprot Proteomes/Eukaryotes/'

for filename in os.listdir(dstDir):
    if filename.endswith(".csv"): 
        if 'with dilution' in filename:
        #if 'without dilution' in filename:
            #print filename
            organismName = filename[:filename.find('(')-1]
            colloquialName = filename[filename.find('(')+1:filename.find(')')]
            
            with open(dstDir + filename, 'rb') as csvfile:
                csvreader = csv.reader(csvfile, delimiter=',')
                headers = next(csvreader, None)[4]
                #print headers
                numOfProts = int(headers[headers.find('Out of ')+len('Out of ')+1:headers.find(' proteins in total')])
                #print numOfProts
                i = 1
                peptides = []
                for row in csvreader:
                    # process each row
                    peptide = row[0]
                    if Counter(peptide).most_common(1)[0][1] >= 8:
                        # ignore almost-SAARs
                        continue
                    if 'X' not in peptide:
                        peptides.append(peptide)
                        i += 1
                    if i > 40: # if i > 40:
                        break
            organisms[organismName + ' (' + colloquialName + ')'] = peptides
            #print colloquialName
            #print organisms[colloquialName]
 
from scipy.spatial import distance

all_peptides = set()
for org in organisms:
    for peptide in organisms[org]:
        all_peptides.add(peptide)

organisms_vectors = {}
for org in organisms:
    organisms_vectors[org] = [1 if peptide in organisms[org] else 0 for peptide in all_peptides]


organisms_order =  ['Arabidopsis Thaliana (arabidopsis)',
                    'Oryza Sativa (rice)',
                    'Saccharomyces Cerevisiae (yeast)',
                    'Candida Glabrata (haploid yeast)',
                    'Caenorhabditis Elegans (worm)',
                    'Biomphalaria Glabrata (mollusca)',
                    'Ciona Intestinalis (sea squirt)',
                    'Drosophila Melanogaster (fly)',
                    'Anopheles Gambiae (mosquito)',
                    'Apis Mellifera (bee)',
                    'Tribolium Castaneum (beetle)',
                    'Danio Rerio (zebrafish)',
                    'Takifugu Rubripes (fugu)',
                    'Tetraodon Nigroviridis (pufferfish)',
                    'Gasterosteus Aculeatus (stickleback)',
                    'Oryzias Latipes (Japanese Medaka)',
                    'Anolis Carolinensis (lizard)',
                    'Xenopus Tropicalis (frog)',
                    'Gallus Gallus (chicken)',
                    'Ornithorhynchus Anatinus (platypus)',
                    'Bos Taurus (cow)',
                    'Canis Familiaris (dog)',
                    'Rattus Norvegicus (rat)',
                    'Macaca Mulatta (rhesus monkey)',
                    'Monodelphis Domestica (opossum)',
                    'Mus Musculus (mouse)',
                    'Pan Troglodytes (chimpanzee)',
                    'Homo Sapiens (human)']
organisms_order = list(reversed(organisms_order))
                    
colloquial_organisms_order =   ['Arabidopsis',
                                'Rice',
                                'S. Cerevisiae',
                                'Candida',
                                'Worm',
                                'Mollusca',
                                'Sea Squirt',
                                'Fly',
                                'Mosquito',
                                'Bee',
                                'Beetle',
                                'Zebrafish',
                                'Fugu',
                                'Pufferfish',
                                'Stickleback',
                                'Japanese Medaka',
                                'Lizard',
                                'Frog',
                                'Chicken',
                                'Platypus',
                                'Cow',
                                'Dog',
                                'Rat',
                                'Opossum',
                                'Mouse',
                                'Rhesus Macaque',
                                'Chimpanzee',
                                'Human']
colloquial_organisms_order = list(reversed(colloquial_organisms_order))
 
dist_matrix = [[distance.hamming(organisms_vectors[org1],organisms_vectors[org2])*len(all_peptides) for org1 in organisms_order ] for org2 in organisms_order]
#dist_matrix = [[distance.hamming(organisms_vectors[org1],organisms_vectors[org2]) for org1 in organisms_order ] for org2 in organisms_order]
#dist_matrix = [[distance.euclidean(organisms_vectors[org1],organisms_vectors[org2]) for org1 in organisms_order ] for org2 in organisms_order]

import matplotlib.pyplot as plt
import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()

uniform_data = np.array(dist_matrix)
mask = np.zeros_like(uniform_data)
mask[np.triu_indices_from(mask)] = True
with sns.axes_style("white"):
    fig, ax = plt.subplots(figsize=(15,15))
    ax = sns.heatmap(uniform_data, xticklabels = colloquial_organisms_order, yticklabels = colloquial_organisms_order, cmap="Reds_r", mask=mask, square = True)
    ax.set_yticklabels(ax.get_yticklabels(), rotation = 0, fontsize = 16)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, fontsize = 16)
    ax.figure.tight_layout()
    ax.figure.subplots_adjust(bottom = 0.25)
    plt.show()
    #plt.savefig('outputs/{}/graph - edges for hamming distance {}.png')
    
"""
sns.set(font_scale=1.8)
uniform_data = np.array(dist_matrix)
mask = np.zeros_like(uniform_data)
mask[np.triu_indices_from(mask)] = True
with sns.axes_style("white"):
    ax = sns.clustermap(uniform_data, xticklabels = colloquial_organisms_order, yticklabels = colloquial_organisms_order, cmap="Reds_r", figsize = (15,15))
    plt.setp(ax.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(ax.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    #ax.set_yticklabels(ax.get_yticklabels(), rotation = 0, fontsize = 18)
    #ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, fontsize = 18)
    #ax.figure.tight_layout()
    #ax.figure.subplots_adjust(bottom = 0.22, right = 0.68)
    plt.show()
"""
