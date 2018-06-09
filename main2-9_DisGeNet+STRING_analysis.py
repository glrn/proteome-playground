from Bio import SeqIO
import csv
import utils
import ast

DisGeNET_file = 'data/DisGeNET/all_gene_disease_associations.tsv'

# Create mapping of diseaseId to diseaseName


disease_id_mapping = dict()

with open(DisGeNET_file) as csvfile:
    reader = csv.DictReader(csvfile, delimiter = '\t')
    for row in reader:
        disease_id_mapping[row['diseaseId']] = row['diseaseName']

print "A total of %d diseases found in DisGeNET" % len(disease_id_mapping)

#####

with open("data/Homo Sapiens - all genes.txt") as f:
    all_human_genes = f.readlines()
all_human_genes = [x.strip() for x in all_human_genes]
all_human_genes = set(all_human_genes)

stringDB = dict.fromkeys(all_human_genes, [])

disease_gene_mapping = dict() # a set of genes for each disease id
MIN_ASSOCIATION_SCORE = 0.1
MIN_GENES_FOR_DISEASE = 10

with open(DisGeNET_file) as csvfile:
    reader = csv.DictReader(csvfile, delimiter = '\t')
    for row in reader:
        if float(row['score']) >= MIN_ASSOCIATION_SCORE:
            diseaseId = row['diseaseId']
            geneSymbol = row['geneSymbol']
            if diseaseId not in disease_gene_mapping:
                disease_gene_mapping[diseaseId] = set()
            disease_gene_mapping[diseaseId].add(geneSymbol)

print "A total of %d diseases with at least %d associated genes (with score > %.2f" % \
      (len([disease for disease in disease_gene_mapping if len(disease_gene_mapping[disease]) >= MIN_GENES_FOR_DISEASE]),
       MIN_GENES_FOR_DISEASE,
       MIN_ASSOCIATION_SCORE)

STRING_INTERACTIONS_FILE = 'human_genome_interactions_min_score_900.csv'

with open(STRING_INTERACTIONS_FILE, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = '\t')
    reader.next() # skip header
    for row in reader:
        gene = row[0]
        interactions = ast.literal_eval(row[1])
        for g in interactions:
            if g not in stringDB:
                stringDB[g] = []
        stringDB[gene] = interactions

print "finished lodaing the stringDB dictionary"
print "number of total genes: %d" % len(stringDB)
print "number of genes with interactions (score > 0.9): %d" % \
      len([x for x in stringDB.values() if x])
print

disease_interactions = dict()
for diseaseId in disease_gene_mapping:
    if len(disease_gene_mapping[diseaseId]) > MIN_GENES_FOR_DISEASE:
        if diseaseId not in disease_interactions:
            disease_interactions[diseaseId] = set()
        for gene in disease_gene_mapping[diseaseId]:
            if gene in stringDB:
                disease_interactions[diseaseId].update(stringDB[gene])


##################

most_frequent_kmers =  ['AAAAAAAAAA','EEEEEEEEEE','QQQQQQQQQQ','PPPPPPPPPP',
                        'SSSSSSSSSS','GGGGGGGGGG','DEEEEEEEEE','EEEEEEEEED',
                        'SAAAAAAAAA','AAAAAAAAAG','PPPPPPPPPL','PPPPPPPPLP',
                        'EEEEEEEEDE','AAAAAAAAAS','EEEEEEEEEG','LPPPPPPPPP',
                        'PPPPPPLPPP','EEEDEEEEEE','PPPPPPPLPP','EDEEEEEEEE',
                        'EEEEDEEEEE','EEDEEEEEEE','PPPPPPPPPA','EEEEEEEDEE',
                        'EEEEEEDEEE','SGGGGGGGGG','APPPPPPPPP','ASSSSSSSSS',
                        'HHHHHHHHHH','QQQQQQQQQP','AAAAAAAAAV','EEEEEDEEEE',
                        'QVKIWFQNRR','EEEEEEEEEA','HQRIHTGEKP','IHTGEKPYKC',
                        'GGGGGGGGGS','CQGDSGGPLV','QPPPPPPPPP','GGGGGGGSGG',
                        'LLLLLLLLLL','VAAAAAAAAA','PPPPLPPPPP','HRRIHTGEKP',
                        'LQQQQQQQQQ','QGDSGGPLVC','GGGGGGSGGG','PSSSSSSSSS']

human_proteome = 'data/Uniprot Proteomes/Eukaryotes/Homo Sapiens (human).fasta'


genes_containing_freq_kmers = set()

for rec in SeqIO.parse(open(human_proteome), 'fasta'):
    seq = str(rec.seq)
    uniqueIdentifier, entryName, proteinName, organismName, geneName = \
        utils.parse_UniProtKB_header(rec.description)
    for freqkmer in most_frequent_kmers:
        if freqkmer in seq:
            genes_containing_freq_kmers.add(geneName)

print "In total found %d genes containing a frequent kmer" % (len(genes_containing_freq_kmers))

for diseaseId in disease_interactions:
    if len(disease_interactions[diseaseId]) >= MIN_GENES_FOR_DISEASE:
        disease_related_genes = disease_interactions[diseaseId].intersection(genes_containing_freq_kmers)
        if len(disease_related_genes) / float(len(disease_interactions[diseaseId])) > 0.01:
            print "diseaseId = %s:\t%d genes,\t%d of them contain a frequent kmer (%.2f)\t[%s]" % \
                  (diseaseId,
                   len(disease_interactions[diseaseId]),
                   len(disease_related_genes),
                   len(disease_related_genes) / float(len(disease_interactions[diseaseId])),
                   disease_id_mapping[diseaseId])