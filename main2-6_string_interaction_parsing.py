import csv
import numpy as np
import ast

with open("data/Homo Sapiens - all genes.txt") as f:
    all_human_genes = f.readlines()
all_human_genes = [x.strip() for x in all_human_genes]
all_human_genes = set(all_human_genes)

stringDB = dict.fromkeys(all_human_genes, [])

#STRING_INTERACTIONS_FILE = 'human_genome_interactions_min_score_900.csv'
STRING_INTERACTIONS_FILE = 'human_genome_interactions_min_score_400.csv'

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

"""
for gene in stringDB:
    if stringDB[gene]:
        degree = len(stringDB[gene])
        degs_of_neighbors = [len(stringDB[n]) for n in stringDB[gene]]
        avg_neighbors_deg = np.average(degs_of_neighbors)
        median_neighbors_deg = np.median(degs_of_neighbors)
    else:
        degree = 0
        avg_neighbors_deg = 0
        median_neighbors_deg = 0
    print '%s\t%d, %.2f, %.2f' % (gene, degree, avg_neighbors_deg, median_neighbors_deg)
"""

all_degrees = [len(stringDB[gene]) for gene in all_human_genes]

print "Degree of gene in homo sapiens: Average -> %.2f\t Median -> %.2f" % (np.average(all_degrees), np.median(all_degrees))
with open('outputs/Homo Sapiens/frequent k10-mers - 2018-02-20-132439.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')
    reader.next() # skip header
    for row in reader:
        kmer = row[0]
        gene_count = int(row[1])
        genes_that_contain_kmer = ast.literal_eval(row[2])
        if gene_count < 10:
            continue
        print "Degree of %d genes that contain %s: Average -> %.2f\t Median -> %.2f" % \
              (gene_count,
               kmer,
               np.average([len(stringDB[g]) for g in genes_that_contain_kmer]),
               np.median([len(stringDB[g]) for g in genes_that_contain_kmer]))

print
print "================"
print

print "Average of average of degrees of neighbors of human genes: %.2f" % \
        np.average([np.average([len(stringDB[g]) for g in stringDB[gene]] if stringDB[gene] else 0) for gene in all_human_genes])
with open('outputs/Homo Sapiens/frequent k10-mers - 2018-02-20-132439.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')
    reader.next() # skip header
    for row in reader:
        kmer = row[0]
        gene_count = int(row[1])
        genes_that_contain_kmer = ast.literal_eval(row[2])
        if gene_count < 10:
            continue
        print "Average of average of degrees of neighDegree of %d genes that contain %s: %.2f" % \
              (gene_count,
               kmer,
               np.average([np.average(
                   [len(stringDB[g]) for g in stringDB[gene]] if stringDB[
                       gene] else 0) for gene in genes_that_contain_kmer]))