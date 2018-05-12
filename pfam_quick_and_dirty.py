from Bio import SeqIO
# local imports
from progressbar import Progressbar
import params
import utils
import sys
import random

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

print "Reading Uniprot file and generating k-mers list for each protein..."
created_protein_names = set() # prevent creation of two similar protein objects
for rec in SeqIO.parse(open("data/Uniprot-Swiss-prot (Homo Sapiens) without fragments.fasta"), 'fasta'):
    seq = str(rec.seq)
    uniqueIdentifier, entryName, proteinName, organismName, geneName = \
        utils.parse_UniProtKB_header(rec.description)
    if geneName == '':
        geneName = uniqueIdentifier
        #print 'Using uncharacterized gene with identifier %s' % uniqueIdentifier
        #print 'Ignoring unknown gene: %s' % rec.description
        #continue
    if geneName in created_protein_names:
        #print 'Ignoring duplicate gene: %s' % rec.description
        continue
    # create a new Protein object
    created_protein_names.add(geneName)
    protein_seq[geneName] = seq
    all_proteins.append(Protein(geneName, seq, params.K))

print
print "Counted k-mer (k=%d) for %d different genes (proteins)." % (params.K, len(all_proteins))

kmer = 'QQQQQQQQQQ'
random.shuffle(all_proteins)
for prot in all_proteins:
    if kmer in prot.seq:
        print '%s (%d)' % (prot.geneName, prot.seq.find(kmer))