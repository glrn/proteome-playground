from Bio import SeqIO
# local imports
from progressbar import Progressbar
import params
import utils
import sys
import os.path

class Protein:
    """
    Hold information of gene name, amino acids, and all k-mers of amino acids.
    """

    def __init__(self, geneName, sequence, k):
        self.geneName = geneName
        self.seq = sequence
        self.kmers = set([self.seq[i:i+k] for i in xrange(len(self.seq)-k)])

def main(proteome_file, similar_diluted):
    all_proteins = list() # all_proteins[i] = PROTEIN()_OBJECT
    protein_seq = dict() # protein_seq[GENE_NAME] = AMINO_ACID_SEQUENCE

    print "Reading Uniprot file and generating k-mers list for each protein..."
    created_protein_names = set() # prevent creation of two similar protein objects
    duplicate_genes_ignored = 0
    for rec in SeqIO.parse(open(proteome_file), 'fasta'):
        seq = str(rec.seq)
        uniqueIdentifier, entryName, proteinName, organismName, geneName = \
            utils.parse_UniProtKB_header(rec.description)
        if geneName == '':
            geneName = uniqueIdentifier
            #print 'Using uncharacterized gene with identifier %s' % uniqueIdentifier
            #print 'Ignoring unknown gene: %s' % rec.description
            #continue
        if geneName in created_protein_names:
            duplicate_genes_ignored += 1
            #print 'Ignoring duplicate gene: %s' % rec.description
            continue
        # create a new Protein object
        created_protein_names.add(geneName)
        protein_seq[geneName] = seq
        all_proteins.append(Protein(geneName, seq, params.K))

    print
    print "Ignored %d duplicate genes." % duplicate_genes_ignored
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
            if similar_diluted:
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
    filename_without_extension = os.path.splitext(os.path.basename(proteome_file))[0]
    dilution_status = 'with dilution' if similar_diluted else 'without dilution'
    outfile = "{0} - frequent k{1}-mers - {2} - {3}.csv".format(filename_without_extension, params.K, dilution_status, timestamp)
    outfile = os.path.join(os.path.dirname(proteome_file),outfile)
    with open(outfile, "wb") as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            writer.writerow(['k-mer', 'number of proteins', 'percentage','all','(Out of %d proteins in total)' % len(all_proteins)])

            for kmer in most_frequenct_kmers:
                total_proteins = len(kmers_frequency[kmer])
                if total_proteins < 5:
                    break
                percentage = round(float(total_proteins) / len(all_proteins), 6)
                geneList = list(kmers_frequency[kmer])
                #geneList = '\r\n'.join(geneList)
                row = [kmer, total_proteins, percentage, geneList]
                writer.writerow(row)

if __name__ == "__main__":
    # eample:  python main2-1_kmer_counting.py \
    #                   "data/UniProt full Proteomes/uniprot-proteome-full Ciona intestinalis (17309).fasta" \
    #                   --similar-diluted
    # output shall be saved in same directory
    proteome_file = sys.argv[1]
    if len(sys.argv) > 2 and sys.argv[2] == '--similar-diluted':
        similar_diluted = True
    else:
        similar_diluted = False

    main(proteome_file, similar_diluted)
