import csv
import ast
# local imports
import params

# Connected components from graph w/ hamming distance <= 3
cc1 = set(['AAAAAAAAAA', 'PAAAAAAAAA', 'GAAAAAAAAA', 'RAAAAAAAAA', 'AAAAAAAAAG', 'AVAAAAAAAA', 'AAAAAAVAAA', 'SAAAAAAAAA', 'AAAAAAAAAV', 'AAAAAAAAAL', 'AAAAAAAAGG', 'VAAAAAAAAA', 'ASAAAAAAAA', 'AAAAAAAAAT', 'AAAAAAAAAS', 'LAAAAAAAAA', 'AAAAAAAASS', 'AAVAAAAAAA', 'TAAAAAAAAA', 'AAAVAAAAAA'])
cc2 = set(['GSGGGGGGGG', 'GGGGGGGGAG', 'AGGGGGGGGG', 'GGGGGGGGSG', 'SGGGGGGGGG', 'GGGGGSGGGG', 'GGGGGGGGGS', 'GGGGGAGGGG', 'GGGGGGGGGG', 'GGGGGGGSGG', 'GGGGSGGGGG', 'GAGGGGGGGG', 'GGSGGGGGGG', 'GGGGGGGGGA', 'GGGGGGSGGG'])
cc3 = set(['PPPPLPPPPP', 'PLPPPPPPPP', 'PPPPPPPPPL', 'PPPPPPPPLP', 'QPPPPPPPPP', 'PPPLPPPPPP', 'LPPPPPPPPP', 'PPPPPLPPPP', 'PPPPPPPPPG', 'APPPPPPPPP', 'PPPPPPLPPP', 'PPPPPPPLPP', 'PPPPPPPPPS', 'PPPPPPPPAP', 'PPPPPPPPPA', 'PPPPPPPPPP', 'PPLPPPPPPP'])
cc4 = set(['QEEEEEEEEE', 'EEEEEDEEEE', 'EEEEEEGEEE', 'EEEGEEEEEE', 'GEEEEEEEEE', 'EEEEEEEEDD', 'EEEEEEEEDE', 'EDEEEEEEEE', 'EEEEEGEEEE', 'DEEEEEEEEE', 'EEEDEEEEEE', 'EEEEEEDEEE', 'AEEEEEEEEE', 'EEDEEEEEEE', 'EEEEDEEEEE', 'SEEEEEEEEE', 'EEEEEEEDED', 'EEEEEEEDEE', 'EEEEEEEEEE', 'EEEEEEEEED', 'EEEEEEEEEG', 'EEEEEEEEEA'])
connected_component = cc4

class Kmer:
    def __init__(self, kmer, gene_count, genes):
        self.kmer = kmer
        self.gene_count = gene_count
        self.genes = genes

# retrieve list of frequent kmers
kmers = {}
with open(params.FREQUENT_KMERS, 'rb') as kmercsv:
    csvreader = csv.reader(kmercsv, delimiter=',')
    next(csvreader, None) # skip header
    for row in csvreader:
        kmer, gene_count, genes = row
        gene_count = int(gene_count)
        genes = ast.literal_eval(genes)
        kmers[kmer] = Kmer(kmer, gene_count, genes)

# build a set of all genes in the connected component
genes_in_connected_component = set()
for kmer in connected_component:
    genes_in_connected_component.update(kmers[kmer].genes)

# write results to file
with open("outputs/connected component 4.txt", 'wb') as file:
    file.write('In total %d genes in connected component:\n' %
               len(genes_in_connected_component))
    for gene in genes_in_connected_component:
        file.write(gene + '\n')
    file.write('------------------\n')
    for kmer in connected_component:
        file.write('Genes that include %s:\n' % kmer)
        file.write('\n'.join(kmers[kmer].genes))
        file.write('\n------------------\n')