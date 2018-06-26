import os
import csv

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
                    if 'X' not in peptide:
                        peptides.append(peptide)
                        i += 1
                    if i > 40:
                        break
                        
            organisms[colloquialName] = peptides
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


organisms_order =  ['human',
                    'chimpanzee',
                    'mouse',
                    'rat',
                    'cow',
                    'dog',
                    'platypus',
                    'chicken',
                    'lizard',                  
                    'fugu',
                    'pufferfish',
                    'beetle',
                    'Japanese Medaka',
                    'frog',
                    'bee',
                    'mosquito',
                    'fly',
                    'worm',
                    'mollusca',
                    'yeast',
                    'haploid yeast']
 
dist_matrix = [[distance.hamming(organisms_vectors[org1],organisms_vectors[org2])*len(all_peptides) for org1 in organisms_order ] for org2 in organisms_order]
#dist_matrix = [[distance.hamming(organisms_vectors[org1],organisms_vectors[org2]) for org1 in organisms_order ] for org2 in organisms_order]
#dist_matrix = [[distance.euclidean(organisms_vectors[org1],organisms_vectors[org2]) for org1 in organisms_order ] for org2 in organisms_order]

import matplotlib.pyplot as plt
import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()
uniform_data = np.array(dist_matrix)
ax = sns.heatmap(uniform_data, xticklabels = organisms_order, yticklabels = organisms_order)
ax.set_yticklabels(ax.get_yticklabels(), rotation = 0, fontsize = 18)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, fontsize = 18)
ax.figure.tight_layout()
plt.show()