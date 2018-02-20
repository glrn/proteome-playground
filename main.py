from Bio import SeqIO
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


all_proteins = list()

print "Reading Uniprot file and generating k-mers list for each protein..."
created_protein_names = set() # prevent creation of two similar protein objects
for rec in SeqIO.parse(open(params.HUMAN_PROTEOME), 'fasta'):
    seq = str(rec.seq)
    uniqueIdentifier, entryName, proteinName, organismName, geneName = \
        utils.parse_UniProtKB_header(rec.description)
    if geneName == '':
        print 'Ignoring unknown gene: %s' % rec.description
        continue
    if geneName in created_protein_names:
        print 'Ignoring duplicate gene: %s' % rec.description
        continue
    # create a new Protein object
    created_protein_names.add(geneName)
    all_proteins.append(Protein(geneName, seq, params.K))

print
print "Counted k-mer (k=%d) for %d different genes (proteins)." % (params.K, len(all_proteins))

for prot in all_proteins:
    print prot.geneName
    print prot.seq
    print prot.kmers


