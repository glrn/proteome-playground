import csv
import sys

### Parameters ###
MIN_ASSOCIATION_SCORE = 0.0
MIN_GENES_FOR_DISEASE = 1
MIN_ENRICHMENT_SCORE = 0.10
##################

try:
    input_file = sys.argv[1]
except:
    print "Error! gene input file exptected."
    print "Example:\tpython script.py example_genes.txt"
    sys.exit()

print "Loading input file %s..." % input_file
genes = open(input_file,'r').read().splitlines()
print "Loaded %d genes from file" % len(genes)


DisGeNET_file = 'DisGeNET_all_gene_disease_associations.tsv'

# Create mapping of diseaseId to diseaseName

disease_id_mapping = dict()

with open(DisGeNET_file) as csvfile:
    reader = csv.DictReader(csvfile, delimiter = '\t')
    for row in reader:
        disease_id_mapping[row['diseaseId']] = row['diseaseName']

print "A total of %d diseases found in DisGeNET" % len(disease_id_mapping)


disease_gene_mapping = dict() # a set of genes for each disease id


with open(DisGeNET_file) as csvfile:
    reader = csv.DictReader(csvfile, delimiter = '\t')
    for row in reader:
        if float(row['score']) >= MIN_ASSOCIATION_SCORE:
            diseaseId = row['diseaseId']
            geneSymbol = row['geneSymbol']
            if diseaseId not in disease_gene_mapping:
                disease_gene_mapping[diseaseId] = set()
            disease_gene_mapping[diseaseId].add(geneSymbol)

print "A total of %d diseases with at least %d associated genes (with score > %.2f)" % \
      (len([disease for disease in disease_gene_mapping if len(disease_gene_mapping[disease]) >= MIN_GENES_FOR_DISEASE]),
       MIN_GENES_FOR_DISEASE,
       MIN_ASSOCIATION_SCORE)

print "Loaded a list of %d genes" % (len(genes))

for diseaseId in disease_gene_mapping:
    if len(disease_gene_mapping[diseaseId]) >= MIN_GENES_FOR_DISEASE:
        disease_related_genes = disease_gene_mapping[diseaseId].intersection(genes)
        if len(disease_related_genes) / float(len(disease_gene_mapping[diseaseId])) > MIN_ENRICHMENT_SCORE:
            print "diseaseId = %s:\tthe disease group contains %d genes, %d of them are genes from the input list (%.2f)\t[Disease: %s]" % \
                  (diseaseId,
                   len(disease_gene_mapping[diseaseId]),
                   len(disease_related_genes),
                   len(disease_related_genes) / float(len(disease_gene_mapping[diseaseId])),
                   disease_id_mapping[diseaseId])

