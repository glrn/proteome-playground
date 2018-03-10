#OUTPUT_DIR = 'Homo Sapiens'
OUTPUT_DIR = 'Drosophila melanogaster'
#OUTPUT_DIR = 'Mus Musculus'

# For kmer counting
#PROTEOME_FILE = 'data/Uniprot-Swiss-prot (Homo Sapiens) without fragments.fasta'
PROTEOME_FILE = 'data/Uniprot-Swiss-prot (Drosophila melanogaster) without fragments.fasta'
#PROTEOME_FILE = 'data/Uniprot-Swiss-prot (Mus Musculus) without fragments.fasta'
K = 10

# For self-repeating kmer counting
MIN_REPETITIONS_IN_PROTEIN = 2
MIN_DIST_BETWEEN_REPETITIONS = K # this means non-overlapping

# For graph drawing
SHOULD_SAVEFIG = False
FREQUENT_KMERS = 'outputs/frequent k10-mers - 2018-02-20-132439.csv'
MIN_HAM_DIST = 3
