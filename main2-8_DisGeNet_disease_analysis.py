from Bio import SeqIO
import csv
import utils

DisGeNET_file = 'data/DisGeNET/all_gene_disease_associations.tsv'

# Create mapping of diseaseId to diseaseName


disease_id_mapping = dict()

with open(DisGeNET_file) as csvfile:
    reader = csv.DictReader(csvfile, delimiter = '\t')
    for row in reader:
        disease_id_mapping[row['diseaseId']] = row['diseaseName']

print "A total of %d diseases found in DisGeNET" % len(disease_id_mapping)



disease_gene_mapping = dict() # a set of genes for each disease id
MIN_ASSOCIATION_SCORE = 0.0
MIN_GENES_FOR_DISEASE = 10
#MIN_GENES_FOR_DISEASE = 1

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

########### new
"""
for i in xrange(50):
    print "================================"
    human_genes = [g.strip() for g in open('data/Homo Sapiens - all genes.txt','r').readlines()]
    import random
    genes_containing_freq_kmers = random.sample(set(human_genes), 100)
    print genes_containing_freq_kmers

    for diseaseId in disease_gene_mapping:
        if len(disease_gene_mapping[diseaseId]) >= MIN_GENES_FOR_DISEASE:
            disease_related_genes = disease_gene_mapping[diseaseId].intersection(genes_containing_freq_kmers)
            if len(disease_related_genes) / float(len(disease_gene_mapping[diseaseId])) > 0.20:
                print "diseaseId = %s:\t%d genes,\t%d of them contain a frequent kmer (%.2f)\t[%s]" % \
                      (diseaseId,
                       len(disease_gene_mapping[diseaseId]),
                       len(disease_related_genes),
                       len(disease_related_genes) / float(len(disease_gene_mapping[diseaseId])),
                       disease_id_mapping[diseaseId])
"""
########### new end


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
#most_frequent_kmers =  ['QQQQQQQQQQ']

human_proteome = 'data/Uniprot Proteomes/Eukaryotes/Homo Sapiens (human).fasta'
"""
freqkmer_gene_mapping = dict()

for rec in SeqIO.parse(open(human_proteome), 'fasta'):
    seq = str(rec.seq)
    uniqueIdentifier, entryName, proteinName, organismName, geneName = \
        utils.parse_UniProtKB_header(rec.description)
    for freqkmer in most_frequent_kmers:
        if freqkmer in seq:
            if freqkmer not in freqkmer_gene_mapping:
                freqkmer_gene_mapping[freqkmer] = set()
            freqkmer_gene_mapping[freqkmer].add(geneName)

for freqkmer in freqkmer_gene_mapping:
    for diseaseId in disease_gene_mapping:
        if len(disease_gene_mapping[diseaseId]) >= MIN_GENES_FOR_DISEASE:
            disease_related_genes = disease_gene_mapping[diseaseId].intersection(freqkmer_gene_mapping[freqkmer])
            if len(disease_related_genes) > 1:
                print "diseaseId = %s:\t%d genes,\t%d of them contain %s" % \
                      (diseaseId,
                       len(disease_gene_mapping[diseaseId]),
                       len(disease_related_genes),
                       freqkmer)
"""

genes_containing_freq_kmers = set()

for rec in SeqIO.parse(open(human_proteome), 'fasta'):
    seq = str(rec.seq)
    uniqueIdentifier, entryName, proteinName, organismName, geneName = \
        utils.parse_UniProtKB_header(rec.description)
    for freqkmer in most_frequent_kmers:
        if freqkmer in seq:
            genes_containing_freq_kmers.add(geneName)

print "In total found %d genes containing a frequent kmer" % (len(genes_containing_freq_kmers))

for diseaseId in disease_gene_mapping:
    if len(disease_gene_mapping[diseaseId]) >= MIN_GENES_FOR_DISEASE:
        disease_related_genes = disease_gene_mapping[diseaseId].intersection(genes_containing_freq_kmers)
        if len(disease_related_genes) / float(len(disease_gene_mapping[diseaseId])) > 0.30:
            print "diseaseId = %s:\t%d genes,\t%d of them contain a frequent kmer (%.2f)\t[%s]" % \
                  (diseaseId,
                   len(disease_gene_mapping[diseaseId]),
                   len(disease_related_genes),
                   len(disease_related_genes) / float(len(disease_gene_mapping[diseaseId])),
                   disease_id_mapping[diseaseId])

# Now make analysis for freqkmer

for freqkmer in most_frequent_kmers:
    print "================================"
    print "===== Frequent kmer is: %s =====" % freqkmer
    genes_containing_freq_kmers = set()

    for rec in SeqIO.parse(open(human_proteome), 'fasta'):
        seq = str(rec.seq)
        uniqueIdentifier, entryName, proteinName, organismName, geneName = \
            utils.parse_UniProtKB_header(rec.description)
        if freqkmer in seq:
            genes_containing_freq_kmers.add(geneName)

    print "In total found %d genes containing a frequent kmer" % (len(genes_containing_freq_kmers))

    for diseaseId in disease_gene_mapping:
        if len(disease_gene_mapping[diseaseId]) >= MIN_GENES_FOR_DISEASE:
            disease_related_genes = disease_gene_mapping[diseaseId].intersection(genes_containing_freq_kmers)
            if len(disease_related_genes) / float(len(disease_gene_mapping[diseaseId])) > 0.20:
                print "diseaseId = %s:\t%d genes,\t%d of them contain a frequent kmer (%.2f)\t[%s]" % \
                      (diseaseId,
                       len(disease_gene_mapping[diseaseId]),
                       len(disease_related_genes),
                       len(disease_related_genes) / float(len(disease_gene_mapping[diseaseId])),
                       disease_id_mapping[diseaseId])
                #print disease_gene_mapping[diseaseId]
                #print disease_related_genes