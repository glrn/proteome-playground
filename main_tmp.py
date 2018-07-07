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


all_proteins = list() # all_proteins[i] = PROTEIN()_OBJECT
protein_seq = dict() # protein_seq[GENE_NAME] = AMINO_ACID_SEQUENCE

print "Reading Uniprot file..."
created_protein_names = set() # prevent creation of two similar protein objects
for rec in SeqIO.parse(open(params.PROTEOME_FILE), 'fasta'):
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

g1 = ['FOXP2', 'EP400', 'SMARCA2', 'THAP11', 'NCOA3', 'MAML2', 'MAML3', 'BMP2K',
     'MED12', 'MED15', 'POU3F2', 'KMT2D', 'MN1', 'IRF2BPL', 'AR', 'KCNN3',
     'FRMPD3', 'NCOA6', 'NFAT5', 'TBP', 'RUNX2', 'MAMLD1', 'HTT', 'ATXN8',
     'ATXN2', 'ATXN1']
g2 = ['ARL6IP4', 'MBTPS2', 'TPRXL', 'TNRC18', 'PPRC1', 'KCNMA1', 'MLLT3',
     'SETD1A', 'ZBTB4', 'ZNF865', 'DACH1', 'TMEM40', 'SRRM2']
g3 = ['ZIC5', 'FNBP4', 'ZFHX4', 'FMNL1', 'SETD1A', 'RAPH1', 'HTT', 'PCLO',
     'ALG13', 'SKOR2', 'LMOD2']
g4 = ['EHMT2', 'ZBTB47', 'CDK11B', 'PELP1', 'CHIC1', 'SCAF1', 'DCAF8L2', 'KAT6B',
     'TTBK1', 'MYT1']
g5 = ['POU3F3', 'SOX21', 'ARX', 'SP8', 'MAZ', 'HOXA13', 'IRF2BPL', 'EVX2',
     'RBM24', 'FOXB2']
g6 = ['POU3F2', 'CAPNS1', 'FUS', 'HOXB3', 'NPAS3', 'AR', 'DACH1', 'ZNF503']
g7 = ['ZNF334', 'ZNF619', 'ZKSCAN5', 'GLI4', 'ZNF473', 'ZSCAN20']
g8 = ['ZNF263', 'ZNF487', 'ZNF268', 'PRDM9', 'ZNF629']
g9 = ['ZNF840P', 'ZSCAN31', 'ZFP62', 'GLI4', 'ZNF561']
g10 = ['ZNF487', 'ZNF629', 'ZNF202', 'PRDM9', 'ZNF577']

all_prots = g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + g9 + g10
genes = set()
genes.update(set(['HMX3', 'B7Z774', 'B4DF88', 'hCG_2039588', 'TLE3', 'BHLHE22', 'Q7KZE5', 'DKFZp762H2012', 'SP5', 'TLE4', 'TBX2', 'A0A1U9X8A2', 'PITX3', 'NKX2-3', 'B7Z4A9', 'B3KQ29', 'B4DP89', 'ZNF703', 'HOXA13', 'PBX3', 'EOMES', 'ZSWIM6', 'B7Z5Q0', 'ZIC2', 'FBXL17', 'SPTLC2', 'FOXB2', 'GSX2', 'SP8', 'Q53EU2', 'PBX1', 'FOXD3', 'IRF2BPL', 'EVX2', 'A8K350', 'Q6NW39', 'B4E2W3', 'GATA6', 'B3KUA2', 'A0A1U9X8A1', 'PBX2', 'NLK', 'B3KM67', 'NKX6-1', 'SOX21', 'RBM47', 'ID4', 'MAZ', 'FLJ20273', 'EGR2', 'RBM24', 'PRDM8', 'ZNF395']))

all_prots = list(genes)
from random import shuffle
shuffle(all_prots)



for p in all_prots:
    print p

# dilute similar proteins
diluted_list = []
pb = Progressbar('Diluting proteins...')
i = 0
for protein_name1 in all_prots:
    i += 1
    pb.update_progress(i, len(all_prots))
    redundantProt = False
    for protein_name in diluted_list:
        try:
            #print '%s: Checking similarity of %s and %s' % (kmer, protein_name, prot.geneName)
            if not utils.proteins_are_dissimilar(protein_name, protein_name1,
                                                 protein_seq[protein_name], protein_seq[protein_name1]):
                redundantProt = True
                break
        except:
            continue
    if not redundantProt:
        diluted_list.append(protein_name1)
    else:
        print 'Ignoring redundant protein: %s' % protein_name1
    redundantProt = False

for p in diluted_list:
    print p

