import csv
import networkx as nx
# local imports
import params

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
    for j in xrange(i, len(G.nodes())):
        import random
        p = random.uniform(0, 1)
        if p > 0.99:
            G.add_edge(G.nodes()[i], G.nodes()[j])

print G.nodes()

# plot graph
import matplotlib.pyplot as plt
color_values = [k.gene_count for k in kmers] # TODO: fix label coloring, there might be a bug here...
nx.draw(G, with_labels=True, node_color=color_values, cmap=plt.cm.Reds_r)
plt.show()
