import csv
import ast
import networkx as nx
# local imports
import params
import utils

class Kmer:
    def __init__(self, kmer, gene_count, genes):
        self.kmer = kmer
        self.gene_count = gene_count
        self.genes = genes

# retrieve list of frequent kmers
kmers = list()
gene_count_dict = {} # gene_count[KMER] = NUMBER_OF_PROTEINS_WITH_THAT_KMER
with open(params.FREQUENT_KMERS, 'rb') as kmercsv:
    csvreader = csv.reader(kmercsv, delimiter=',')
    next(csvreader, None) # skip header
    for row in csvreader:
        kmer, gene_count, genes = row
        gene_count = int(gene_count)
        genes = ast.literal_eval(genes)
        if gene_count >= params.MIN_PROTEINS_FOR_KMER_GRAPH:
            kmers.append(Kmer(kmer, gene_count, genes))
            gene_count_dict[kmer] = gene_count

# draw kmer graph
G = nx.Graph()
G.add_nodes_from([k.kmer for k in kmers])

for i in xrange(len(G.nodes())):
    for j in xrange(i+1, len(G.nodes())):
        kmerA, kmerB = G.nodes()[i], G.nodes()[j]
        # decide if there should be an edge between kmerA and kmerB
        if utils.hamdist(kmerA, kmerB) <= params.MIN_HAM_DIST:
            G.add_edge(kmerA, kmerB)

# output info on connected components
for c in nx.connected_components(G):
    if len(c) > 10:
        print c

# plot graph
import matplotlib.pyplot as plt
color_values = [90 - gene_count_dict[kmer] for kmer in G.nodes()]
plt.figure(figsize=(20, 15))
nx.draw(G, with_labels=True, node_color=color_values, vmin=3, vmax=86,cmap=plt.cm.Reds_r)

if params.SHOULD_SAVEFIG:
    plt.savefig('outputs/{}/graph - edges for hamming distance {}.png'.format(
        params.OUTPUT_DIR,params.MIN_HAM_DIST))
else:
    plt.show()