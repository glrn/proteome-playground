# For kmer counting
HUMAN_PROTEOME = 'data/Uniprot-Swiss-prot (Homo Sapiens) without fragments.fasta'
K = 10

# For self-repeating kmer counting
MIN_REPETITIONS_IN_PROTEIN = 2
MIN_DIST_BETWEEN_REPETITIONS = K # this means non-overlapping

# For graph drawing
SHOULD_SAVEFIG = False
FREQUENT_KMERS = 'outputs/frequent k10-mers - 2018-02-20-132439.csv'
MIN_HAM_DIST = 3
