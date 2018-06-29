
list_of_species =  [#'Arabidopsis Thaliana (arabidopsis)',
                    #'Oryza Sativa (rice)',
                    #'Vitis Vinifera (grape)',
                    'Saccharomyces Cerevisiae (yeast)',
                    'Candida Glabrata (haploid yeast)',
                    'Caenorhabditis Elegans (worm)',
                    'Biomphalaria Glabrata (mollusca)',
                    'Ciona Intestinalis (sea squirt)',
                    'Drosophila Melanogaster (fly)',
                    'Anopheles Gambiae (mosquito)',
                    'Apis Mellifera (bee)',
                    'Tribolium Castaneum (beetle)',
                    'Danio Rerio (zebrafish)',
                    'Takifugu Rubripes (fugu)',
                    'Tetraodon Nigroviridis (pufferfish)',
                    'Gasterosteus Aculeatus (stickleback)',
                    'Oryzias Latipes (Japanese Medaka)',
                    'Anolis Carolinensis (lizard)',
                    'Xenopus Tropicalis (frog)',
                    'Gallus Gallus (chicken)',
                    'Ornithorhynchus Anatinus (platypus)',
                    'Bos Taurus (cow)',
                    'Canis Familiaris (dog)',
                    'Rattus Norvegicus (rat)',
                    #'Macaca Mulatta (rhesus monkey)',
                    'Monodelphis Domestica (opossum)',
                    'Mus Musculus (mouse)',
                    'Pan Troglodytes (chimpanzee)',
                    'Homo Sapiens (human)']

data_path = '../data/Uniprot Proteomes/Eukaryotes/'                    
               
files_without_dilution = {}               
#files_without_dilution['Arabidopsis Thaliana (arabidopsis)']     = 
#files_without_dilution['Oryza Sativa (rice)']                    = 
#files_without_dilution['Vitis Vinifera (grape)']                 = 
files_without_dilution['Saccharomyces Cerevisiae (yeast)']       = data_path + 'Saccharomyces Cerevisiae (yeast) - frequent k10-mers - without dilution - 2018-05-06-212054.csv'
files_without_dilution['Candida Glabrata (haploid yeast)']       = data_path + 'Candida Glabrata (haploid yeast) - frequent k10-mers - without dilution - 2018-05-06-212050.csv'
files_without_dilution['Caenorhabditis Elegans (worm)']          = data_path + 'Caenorhabditis Elegans (worm) - frequent k10-mers - without dilution - 2018-05-06-212302.csv'
files_without_dilution['Biomphalaria Glabrata (mollusca)']       = data_path + 'Biomphalaria Glabrata (mollusca) - frequent k10-mers - without dilution - 2018-05-06-212325.csv'
files_without_dilution['Ciona Intestinalis (sea squirt)']        = data_path + 'Ciona Intestinalis (sea squirt) - frequent k10-mers - without dilution - 2018-05-06-212255.csv'
files_without_dilution['Drosophila Melanogaster (fly)']          = data_path + 'Drosophila Melanogaster (fly) - frequent k10-mers - without dilution - 2018-05-06-212313.csv'
files_without_dilution['Anopheles Gambiae (mosquito)']           = data_path + 'Anopheles Gambiae (mosquito) - frequent k10-mers - without dilution - 2018-05-06-212237.csv'
files_without_dilution['Apis Mellifera (bee)']                   = data_path + 'Apis Mellifera (bee) - frequent k10-mers - without dilution - 2018-05-06-212204.csv'
files_without_dilution['Tribolium Castaneum (beetle)']           = data_path + 'Tribolium Castaneum (beetle) - frequent k10-mers - without dilution - 2018-05-06-212250.csv'
files_without_dilution['Danio Rerio (zebrafish)']                = data_path + 'Danio Rerio (zebrafish) - frequent k10-mers - without dilution - 2018-05-06-212425.csv'
files_without_dilution['Takifugu Rubripes (fugu)']               = data_path + 'Takifugu Rubripes (fugu) - frequent k10-mers - without dilution - 2018-05-06-212414.csv'
files_without_dilution['Tetraodon Nigroviridis (pufferfish)']    = data_path + 'Tetraodon Nigroviridis (pufferfish) - frequent k10-mers - without dilution - 2018-05-06-212436.csv'
files_without_dilution['Gasterosteus Aculeatus (stickleback)']   = data_path + 'Gasterosteus Aculeatus (stickleback) - frequent k10-mers - without dilution - 2018-05-06-212346.csv'
files_without_dilution['Oryzias Latipes (Japanese Medaka)']      = data_path + 'Oryzias Latipes (Japanese Medaka) - frequent k10-mers - without dilution - 2018-05-06-212339.csv'
files_without_dilution['Anolis Carolinensis (lizard)']           = data_path + 'Anolis Carolinensis (lizard) - frequent k10-mers - without dilution - 2018-05-06-212307.csv'
files_without_dilution['Xenopus Tropicalis (frog)']              = data_path + 'Xenopus Tropicalis (frog) - frequent k10-mers - without dilution - 2018-05-06-212330.csv'
files_without_dilution['Gallus Gallus (chicken)']                = data_path + 'Gallus Gallus (chicken) - frequent k10-mers - without dilution - 2018-05-06-212319.csv'
files_without_dilution['Ornithorhynchus Anatinus (platypus)']    = data_path + 'Ornithorhynchus Anatinus (platypus) - frequent k10-mers - without dilution - 2018-05-06-212255.csv'
files_without_dilution['Bos Taurus (cow)']                       = data_path + 'Bos Taurus (cow) - frequent k10-mers - without dilution - 2018-05-06-212344.csv'
files_without_dilution['Canis Familiaris (dog)']                 = data_path + 'Canis Familiaris (dog) - frequent k10-mers - without dilution - 2018-05-06-212342.csv'
files_without_dilution['Rattus Norvegicus (rat)']                = data_path + 'Rattus Norvegicus (rat) - frequent k10-mers - without dilution - 2018-05-06-212322.csv'
#files_without_dilution['Macaca Mulatta (rhesus monkey)']         = 
files_without_dilution['Monodelphis Domestica (opossum)']        = data_path + 'Monodelphis Domestica (opossum) - frequent k10-mers - without dilution - 2018-05-06-212337.csv'
files_without_dilution['Mus Musculus (mouse)']                   = data_path + 'Mus Musculus (mouse) - frequent k10-mers - without dilution - 2018-05-06-212349.csv'
files_without_dilution['Pan Troglodytes (chimpanzee)']           = data_path + 'Pan Troglodytes (chimpanzee) - frequent k10-mers - without dilution - 2018-05-06-212353.csv'
files_without_dilution['Homo Sapiens (human)']                   = data_path + 'Homo Sapiens (human) - frequent k10-mers - without dilution - 2018-05-06-212428.csv'

files_with_dilution = {}
#files_with_dilution['Arabidopsis Thaliana (arabidopsis)']        = 
#files_with_dilution['Oryza Sativa (rice)']                       = 
#files_with_dilution['Vitis Vinifera (grape)']                    = 
files_with_dilution['Saccharomyces Cerevisiae (yeast)']          = data_path + 'Saccharomyces Cerevisiae (yeast) - frequent k10-mers - with dilution - 2018-05-06-232625.csv'
files_with_dilution['Candida Glabrata (haploid yeast)']          = data_path + 'Candida Glabrata (haploid yeast) - frequent k10-mers - with dilution - 2018-05-06-212707.csv'
files_with_dilution['Caenorhabditis Elegans (worm)']             = data_path + 'Caenorhabditis Elegans (worm) - frequent k10-mers - with dilution - 2018-05-06-214145.csv'
files_with_dilution['Biomphalaria Glabrata (mollusca)']          = data_path + 'Biomphalaria Glabrata (mollusca) - frequent k10-mers - with dilution - 2018-05-06-215927.csv'
files_with_dilution['Ciona Intestinalis (sea squirt)']           = data_path + 'Ciona Intestinalis (sea squirt) - frequent k10-mers - with dilution - 2018-06-29-013933.csv'
files_with_dilution['Drosophila Melanogaster (fly)']             = data_path + 'Drosophila Melanogaster (fly) - frequent k10-mers - with dilution - 2018-05-07-134053.csv'
files_with_dilution['Anopheles Gambiae (mosquito)']              = data_path + 'Anopheles Gambiae (mosquito) - frequent k10-mers - with dilution - 2018-05-07-055151.csv'
files_with_dilution['Apis Mellifera (bee)']                      = data_path + 'Apis Mellifera (bee) - frequent k10-mers - with dilution - 2018-05-06-222607.csv'
files_with_dilution['Tribolium Castaneum (beetle)']              = data_path + 'Tribolium Castaneum (beetle) - frequent k10-mers - with dilution - 2018-05-07-001436.csv'
files_with_dilution['Danio Rerio (zebrafish)']                   = data_path + 'Danio Rerio (zebrafish) - frequent k10-mers - with dilution - 2018-06-29-050742.csv'
files_with_dilution['Takifugu Rubripes (fugu)']                  = data_path + 'Takifugu Rubripes (fugu) - frequent k10-mers - with dilution - 2018-05-07-200816.csv'
files_with_dilution['Tetraodon Nigroviridis (pufferfish)']       = data_path + 'Tetraodon Nigroviridis (pufferfish) - frequent k10-mers - with dilution - 2018-05-09-063318.csv'
files_with_dilution['Gasterosteus Aculeatus (stickleback)']      = data_path + 'Gasterosteus Aculeatus (stickleback) - frequent k10-mers - with dilution - 2018-06-29-053638.csv'
files_with_dilution['Oryzias Latipes (Japanese Medaka)']         = data_path + 'Oryzias Latipes (Japanese Medaka) - frequent k10-mers - with dilution - 2018-05-07-005525.csv'
files_with_dilution['Anolis Carolinensis (lizard)']              = data_path + 'Anolis Carolinensis (lizard) - frequent k10-mers - with dilution - 2018-05-06-233027.csv'
files_with_dilution['Xenopus Tropicalis (frog)']                 = data_path + 'Xenopus Tropicalis (frog) - frequent k10-mers - with dilution - 2018-05-07-063004.csv'
files_with_dilution['Gallus Gallus (chicken)']                   = data_path + 'Gallus Gallus (chicken) - frequent k10-mers - with dilution - 2018-05-07-021300.csv'
files_with_dilution['Ornithorhynchus Anatinus (platypus)']       = data_path + 'Ornithorhynchus Anatinus (platypus) - frequent k10-mers - with dilution - 2018-05-07-005547.csv'
files_with_dilution['Bos Taurus (cow)']                          = data_path + 'Bos Taurus (cow) - frequent k10-mers - with dilution - 2018-05-07-003323.csv'
files_with_dilution['Canis Familiaris (dog)']                    = data_path + 'Canis Familiaris (dog) - frequent k10-mers - with dilution - 2018-05-07-003054.csv'
files_with_dilution['Rattus Norvegicus (rat)']                   = data_path + 'Rattus Norvegicus (rat) - frequent k10-mers - with dilution - 2018-05-07-050120.csv'
#files_with_dilution['Macaca Mulatta (rhesus monkey)']            = 
files_with_dilution['Monodelphis Domestica (opossum)']           = data_path + 'Monodelphis Domestica (opossum) - frequent k10-mers - with dilution - 2018-06-29-055832.csv'
files_with_dilution['Mus Musculus (mouse)']                      = data_path + 'Mus Musculus (mouse) - frequent k10-mers - with dilution - 2018-05-07-080946.csv'
files_with_dilution['Pan Troglodytes (chimpanzee)']              = data_path + 'Pan Troglodytes (chimpanzee) - frequent k10-mers - with dilution - 2018-05-07-234939.csv'
files_with_dilution['Homo Sapiens (human)']                      = data_path + 'Homo Sapiens (human) - frequent k10-mers - with dilution - 2018-05-09-145033.csv'


import csv
from collections import Counter

peptides = {}

for organism in list_of_species:
    peptidesWithoutDilution = []
    peptidesWithDilution = []
    peptidesWithDilutionWithoutSAARs = []
    
    with open(files_without_dilution[organism], 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        headers = next(csvreader, None)[4]
        #print headers
        i = 1
        for row in csvreader:
            # process each row
            peptide = row[0]
            numOfProts = int(row[1])
            if 'X' not in peptide:
                peptidesWithoutDilution.append((peptide,numOfProts))
                i += 1
            if i > 40:
                break
                
    with open(files_with_dilution[organism], 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        headers = next(csvreader, None)[4]
        #print headers
        i = 1
        for row in csvreader:
            # process each row
            peptide = row[0]
            numOfProts = int(row[1])
            #if Counter(peptide).most_common(1)[0][1] >= 8:
                # ignore almost-SAARs
            #    continue
            if 'X' not in peptide:
                peptidesWithDilution.append((peptide,numOfProts))
                i += 1
            if i > 40:
                break
                
    with open(files_with_dilution[organism], 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        headers = next(csvreader, None)[4]
        #print headers
        i = 1
        for row in csvreader:
            # process each row
            peptide = row[0]
            numOfProts = int(row[1])
            if Counter(peptide).most_common(1)[0][1] >= 8:
                # ignore almost-SAARs
                continue
            if 'X' not in peptide:
                peptidesWithDilutionWithoutSAARs.append((peptide,numOfProts))
                i += 1
            if i > 40:
                break
                
    peptides[organism] = (peptidesWithoutDilution, peptidesWithDilution, peptidesWithDilutionWithoutSAARs)
    
for org in list_of_species:
    print org
    nonDiluted =            [p[0] for p in peptides[org][0]]
    diluted =               [p[0] for p in peptides[org][1]]
    dilutedWithoutSAARs =   [p[0] for p in peptides[org][2]]
    #print nonDiluted
    #print diluted
    #print set(diluted).intersection(set(nonDiluted))
    print len(set(diluted).intersection(set(nonDiluted)))
    print len(set(diluted).intersection(set(dilutedWithoutSAARs)))
    print
