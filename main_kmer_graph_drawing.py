import csv
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
with open(params.FREQUENT_KMERS, 'rb') as kmercsv:
    csvreader = csv.reader(kmercsv, delimiter=',')
    next(csvreader, None) # skip header
    for row in csvreader:
        kmer, gene_count, genes = row
        gene_count = int(gene_count)
        if gene_count >= 10:
            kmers.append(Kmer(kmer, int(gene_count), genes))

# draw kmer graph
G = nx.Graph()
G.add_nodes_from([k.kmer for k in kmers])

for i in xrange(len(G.nodes())):
    for j in xrange(i+1, len(G.nodes())):
        kmerA, kmerB = G.nodes()[i], G.nodes()[j]
        # decide if there should be an edge between kmerA and kmerB
        if utils.hamdist(kmerA, kmerB) <=5:
            G.add_edge(kmerA, kmerB)

print G.nodes()

# plot graph
import matplotlib.pyplot as plt
color_values = [k.gene_count for k in kmers] # TODO: fix label coloring, there might be a bug here...
nx.draw(G, with_labels=True, node_color=color_values, cmap=plt.cm.Reds_r)
plt.show()
