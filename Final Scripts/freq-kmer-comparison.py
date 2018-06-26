import csv
from scipy.spatial import distance

lines = open('Frequent k-mers comparison.csv', 'rb').readlines()
lines = [line.rstrip().split(',') for line in lines]

organisms = {}

# read organism name
for i in xrange(len(lines[1])):
    name = lines[1][i]
    organisms[name] = [lines[j+2][i] for j in xrange(len(lines)-2)]
    
#for org in ['Homo sapiens','Pan troglodytes','Mus musculus','Caenorhabditis elegans']:
#    print org
#    for p in organisms[org]:
#        print p
#    print
    
all_peptides = set()
for org in organisms:
    for peptide in organisms[org]:
        all_peptides.add(peptide)
        
#print len(all_peptides)

organisms_vectors = {}
for org in organisms:
    organisms_vectors[org] = [1 if peptide in organisms[org] else 0 for peptide in all_peptides]
    
for org in ['Human','chimpanzee']:
    print org
    print organisms_vectors[org]
    
#print distance.hamming(organisms_vectors['Homo sapiens'],organisms_vectors['Pan troglodytes'])
#print distance.hamming(organisms_vectors['Homo sapiens'],organisms_vectors['Caenorhabditis elegans'])

organisms_order = ['Human',
                   'chimpanzee',
                   'house mouse',
                   'brown rat',
                   'cow',
                   'opossum',
                   'platypus',
                   'chicken',
                   'zebra finch',
                   'lizard',
                   'zebrafish',
                   'bee',
                   'D. melanogaster',
                   'D. yakuba',
                   'sea urchin',
                   'ciona',
                   'worm',
                   'C. elegans',
                   'sea anemone',
                   'yeast',
                   'Asian rice']

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



