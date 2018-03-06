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
        all_kmers = list([self.seq[i:i+k] for i in xrange(len(self.seq)-k)])
        self.kmers = {} # self.kmers is a dictionary that counts number of
                        # occurrences for each kmer in the sequence
        for kmer in all_kmers:
            if kmer not in self.kmers:
                self.kmers[kmer] = 0
            self.kmers[kmer] += 1


all_proteins = list() # all_proteins[i] = PROTEIN()_OBJECT
protein_seq = dict() # protein_seq[GENE_NAME] = AMINO_ACID_SEQUENCE

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
    protein_seq[geneName] = seq
    all_proteins.append(Protein(geneName, seq, params.K))

print
print "Counted k-mer (k=%d) for %d different genes (proteins)." % (params.K, len(all_proteins))

pb = Progressbar('Generating frequency dictionary for k-mers that appear more than once in protein')
kmers_frequency = dict() # track popularity of kmer accross all proteins
i = 0
for prot in all_proteins:
    i += 1
    pb.update_progress(i, len(all_proteins))

    # consider only kmers that appear more than once in their protein,
    self_repeating_kmers = \
        [kmer for kmer in prot.kmers if prot.kmers[kmer] >= params.MIN_REPETITIONS_IN_PROTEIN]

    # remove SAARs
    # self_repeating_kmers = [kmer for kmer in self_repeating_kmers if len(set([c for c in kmer])) > 1]

    # consider only kmers where there are two occurrences that are distanced
    # from each other (in particular, non overlapping)
    self_repeating_kmers  = [kmer for kmer in self_repeating_kmers if
                             prot.seq.find(kmer) + params.MIN_DIST_BETWEEN_REPETITIONS <
                                prot.seq.rfind(kmer)]

    for kmer in self_repeating_kmers:

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
outfile = 'outputs/frequent self-repeating k{}-mers - {}.csv'.format(params.K, timestamp)
with open(outfile, "wb") as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(['k-mer','number of proteins','all'])

        for kmer in most_frequenct_kmers:
            total_proteins = len(kmers_frequency[kmer])
            if total_proteins < 5:
                break
            geneList = list(kmers_frequency[kmer])
            #geneList = '\r\n'.join(geneList)
            row = [kmer, total_proteins, geneList]
            writer.writerow(row)

