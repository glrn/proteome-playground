from Bio import SeqIO
# local imports
from progressbar import Progressbar
import params
import utils

class Protein:
    """
    Hold information of gene name, amino acids, and all k-mers of amino acids.
    """

    def __init__(self, geneName, sequence, k):
        self.geneName = geneName
        self.seq = sequence
        self.kmers = set([self.seq[i:i+k] for i in xrange(len(self.seq)-k)])


all_proteins = list() # all_proteins[i] = PROTEIN()_OBJECT
protein_seq = dict() # protein_seq[GENE_NAME] = AMINO_ACID_SEQUENCE

#PROTEOME_FILE = 'data/NCBI-RefSeq Bos taurus (13192 proteins).fasta'
#OUTPUT_DIR = 'NCBI/Bos taurus' # Done
#PROTEOME_FILE = 'data/NCBI-RefSeq Canis lupus familiaris (58484 proteins).fasta'
#OUTPUT_DIR = 'NCBI/Canis lupus'
#PROTEOME_FILE = 'data/NCBI-RefSeq Danio rerio (52743 proteins).fasta'
#OUTPUT_DIR = 'NCBI/Danio rerio'
#PROTEOME_FILE = 'data/NCBI-RefSeq Gallus gallus (6396 proteins).fasta'
#OUTPUT_DIR = 'NCBI/Gallus gallus' # Done
PROTEOME_FILE = 'data/NCBI-RefSeq Mus musculus (33367 proteins).fasta'
OUTPUT_DIR = 'NCBI/Mus musculus'
#PROTEOME_FILE = 'data/NCBI-RefSeq Pan troglodytes (80607 proteins).fasta'
#OUTPUT_DIR = 'NCBI/Pan troglodytes'
#PROTEOME_FILE = 'data/NCBI-RefSeq Rattus Norvegicus (17438 proteins).fasta'
#OUTPUT_DIR = 'NCBI/Rattus Norvegicus' # Done
#PROTEOME_FILE = 'data/NCBI-RefSeq Xenopus laevis (10950 proteins).fasta'
#OUTPUT_DIR = 'NCBI/Xenopus laevis' # Done
#PROTEOME_FILE = 'data/NCBI-RefSeq Xenopus tropicalis (8579 proteins).fasta'
#OUTPUT_DIR = 'NCBI/Xenopus tropicalis' # Done

print "Reading Uniprot file and generating k-mers list for each protein..."
created_protein_names = set() # prevent creation of two similar protein objects
for rec in SeqIO.parse(open(PROTEOME_FILE), 'fasta'):
    seq = str(rec.seq)
    organismName, geneName = utils.parse_NCBIRefSeq_header(rec.description)
    if geneName == '':
        print 'Ignoring unknown gene: %s' % rec.description
        continue
    if geneName in created_protein_names:
        print 'Ignoring duplicate gene: %s' % rec.description
        continue
    # create a new Protein object
    created_protein_names.add(geneName)
    protein_seq[geneName] = seq
    all_proteins.append(Protein(geneName, seq, params.K))

print
print "Counted k-mer (k=%d) for %d different genes (proteins)." % (params.K, len(all_proteins))

pb = Progressbar('Generating frequency dictionary for k-mers')
kmers_frequency = dict() # track popularity of kmer accross all proteins
i = 0
for prot in all_proteins:
    i += 1
    pb.update_progress(i, len(all_proteins))

    for kmer in prot.kmers:
        if kmer not in kmers_frequency:
            kmers_frequency[kmer] = set()

        # add the new protein only if it's dissimilar enough from all other
        # proteins that were already added and contain this kmer
        protein_names = kmers_frequency[kmer] # list of all prots that share this kmer
        redundantProt = False
        for protein_name in protein_names:
            #print '%s: Checking similarity of %s and %s' % (kmer, protein_name, prot.geneName)
            if not utils.proteins_are_dissimilar(protein_name, prot.geneName,
                                                 protein_seq[protein_name], prot.seq):
                redundantProt = True
                break
        if not redundantProt:
            kmers_frequency[kmer].add(prot.geneName)
        redundantProt = False

print "Sorting frequent k-mers by frequency..."
most_frequenct_kmers = sorted(kmers_frequency, key=lambda k: len(kmers_frequency[k]), reverse=True)

print "Writing results to file..."
import datetime, time, csv
timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H%M%S')
outfile = 'outputs/{}/frequent k{}-mers - {}.csv'.format(OUTPUT_DIR, params.K, timestamp)
with open(outfile, "wb") as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(['k-mer','number of proteins','all', '(Out of %d proteins in total)' % len(all_proteins)])

        for kmer in most_frequenct_kmers:
            total_proteins = len(kmers_frequency[kmer])
            if total_proteins < 5:
                break
            geneList = list(kmers_frequency[kmer])
            #geneList = '\r\n'.join(geneList)
            row = [kmer, total_proteins, geneList]
            writer.writerow(row)

