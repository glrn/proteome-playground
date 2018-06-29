import os
import csv
from collections import Counter

organisms = {}

dstDir = '../data/Uniprot Proteomes/Eukaryotes/'

for filename in os.listdir(dstDir):
    if filename.endswith(".csv"): 
        #if 'with dilution' in filename:
        if 'without dilution' in filename:
            #print filename
            organismName = filename[:filename.find('(')-1]
            colloquialName = filename[filename.find('(')+1:filename.find(')')]
            fullOrgName = '%s (%s)' % (organismName, colloquialName)
            
            with open(dstDir + filename, 'rb') as csvfile:
                csvreader = csv.reader(csvfile, delimiter=',')
                headers = next(csvreader, None)[4]
                #print headers
                numOfProts = int(headers[headers.find('Out of ')+len('Out of '):headers.find(' proteins in total')])
                #print numOfProts
                i = 1
                peptides = []
                for row in csvreader:
                    # process each row
                    peptide = row[0]
                    numProt = row[1]
                    partOfProts = row[2]
                    #if Counter(peptide).most_common(1)[0][1] >= 8:
                        # ignore almost-SAARs
                    #    continue
                    if 'X' not in peptide:
                        peptides.append((peptide,numProt,partOfProts))
                        i += 1
                    if i > 40:
                        break
                        
            organisms[colloquialName] = peptides
            print fullOrgName
            print numOfProts
            print organisms[colloquialName]
            print
            print
 
