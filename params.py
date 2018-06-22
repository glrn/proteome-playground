#Organism name
ORGANISM = 'Homo Sapiens'
#ORGANISM = 'Mus Musculus'
#ORGANISM = 'Drosophila Melanogaster'
#ORGANISM = 'Danio Rerio'
#ORGANISM = 'Caenorhabditis Elegans'
#ORGANISM = 'Rattus Norvegicus'
#ORGANISM = 'Bos Taurus'

OUTPUT_DIR = 'outputs/' + ORGANISM

# For kmer counting
PROTEOME_FILE = 'data/Uniprot-Swiss-prot ({}) without fragments.fasta'.format(ORGANISM)
K = 10

# For self-repeating kmer counting
MIN_REPETITIONS_IN_PROTEIN = 2
MIN_DIST_BETWEEN_REPETITIONS = K # this means non-overlapping

# For graph drawing
SHOULD_SAVEFIG = False # True
MIN_HAM_DIST = 3

if ORGANISM == 'Homo Sapiens':
    FREQUENT_KMERS = 'outputs/Homo Sapiens/frequent k10-mers - 2018-02-20-132439.csv'
    MIN_PROTEINS_FOR_KMER_GRAPH = 10
elif ORGANISM == 'Mus Musculus':
    FREQUENT_KMERS = 'outputs/Mus Musculus/frequent k10-mers - 2018-03-10-155410.csv'
    MIN_PROTEINS_FOR_KMER_GRAPH = 8
elif ORGANISM == 'Drosophila Melanogaster':
    FREQUENT_KMERS = 'outputs/Drosophila Melanogaster/frequent k10-mers - 2018-03-10-151806.csv'
    MIN_PROTEINS_FOR_KMER_GRAPH = 6
elif ORGANISM == 'Danio Rerio':
    FREQUENT_KMERS = 'outputs/Danio Rerio/frequent k10-mers - 2018-03-10-160131.csv'
    MIN_PROTEINS_FOR_KMER_GRAPH = 4
elif ORGANISM == 'Caenorhabditis Elegans':
    FREQUENT_KMERS = 'outputs/Caenorhabditis Elegans/frequent k10-mers - 2018-03-10-155559.csv'
    MIN_PROTEINS_FOR_KMER_GRAPH = 4

