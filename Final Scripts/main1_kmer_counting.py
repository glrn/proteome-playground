from Bio import SeqIO
# local imports
from progressbar import Progressbar
import params
import utils
import sys

class Protein:
    """
    Hold information of gene name, amino acids, and all k-mers of amino acids.
    """

    def __init__(self, geneName, sequence, k):
        self.geneName = geneName
        self.seq = sequence
        self.kmers = set([self.seq[i:i+k] for i in xrange(len(self.seq)-k)])

def main(proteome_file, output_dir):
    all_proteins = list() # all_proteins[i] = PROTEIN()_OBJECT
    protein_seq = dict() # protein_seq[GENE_NAME] = AMINO_ACID_SEQUENCE

    print "Reading Uniprot file and generating k-mers list for each protein..."
    created_protein_names = set() # prevent creation of two similar protein objects
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
            print 'Ignoring duplicate gene: %s' % rec.description
            continue
        # create a new Protein object
        created_protein_names.add(geneName)
        protein_seq[geneName] = seq
        all_proteins.append(Protein(geneName, seq, params.K))

    print
    print "Counted k-mer (k=%d) for %d different genes (proteins)." % (params.K, len(all_proteins))

    pb = Progressbar('Generating frequency dictionary for k-mers')
    skipped_prots = 0
    kmers_frequency = dict() # track popularity of kmer accross all proteins
    i = 0
    for prot in all_proteins:
        i += 1
        pb.update_progress(i, len(all_proteins))

        """
        if prot.geneName.startswith('ZNF') or prot.geneName.startswith('ZF'):
            skipped_prots += 1
            continue
        if prot.geneName.startswith('OR'):
            skipped_prots += 1
            continue
        if prot.geneName.startswith('HOX'):
            skipped_prots += 1
            continue
        if prot.geneName.startswith('IGKV'):
            skipped_prots += 1
            continue
        """

        for kmer in prot.kmers:
            if kmer not in kmers_frequency:
                kmers_frequency[kmer] = set()
            kmers_frequency[kmer].add(prot.geneName)

    print "Sorting frequent k-mers by frequency..."
    most_frequenct_kmers = sorted(kmers_frequency, key=lambda k: len(kmers_frequency[k]), reverse=True)
    import datetime, time, csv
    timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H%M%S')
    outfile = '{}/frequent k{}-mers - {}.csv'.format(output_dir, params.K, timestamp)
    with open(outfile, "wb") as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            writer.writerow(['k-mer', 'number of proteins', 'percentage of total', 'all','(Out of %d proteins in total)' % (len(all_proteins) - skipped_prots)])

            for kmer in most_frequenct_kmers:
                total_proteins = len(kmers_frequency[kmer])
                if total_proteins < 10:
                    break
                percentage = "{0:.4f}".format(float(total_proteins)/(len(all_proteins) - skipped_prots))
                geneList = list(kmers_frequency[kmer])
                #geneList = '\r\n'.join(geneList)
                row = [kmer, total_proteins, percentage, geneList]
                writer.writerow(row)

if __name__ == "__main__":
    proteome_file = '../data/Uniprot-Swiss-prot (Homo Sapiens) without fragments.fasta'
    output_dir = '.'
    main(proteome_file, output_dir)
