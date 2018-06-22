import csv, ast
from itertools import chain

CSV_FILE = 'frequent k10-mers - 2018-06-16-092428.csv'

genes_that_share_peptide_with_at_least_x_others = dict()

with open(CSV_FILE, 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    next(spamreader) # skip header
    for row in spamreader:
        peptide = row[0]
        genes = ast.literal_eval(row[3])
        num_genes = len(genes)
        #if num_genes < 5:
        #    continue
        if num_genes not in genes_that_share_peptide_with_at_least_x_others:
            genes_that_share_peptide_with_at_least_x_others[num_genes] = set()
        genes_that_share_peptide_with_at_least_x_others[num_genes].update(genes)

for i in range(1, 183+1):
    print '%d\t%d' % (i, len(set(chain.from_iterable([genes_that_share_peptide_with_at_least_x_others.get(j, []) for j in range(i, 185+1)]))))