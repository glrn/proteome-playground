from Bio import SeqIO
import csv
import utils

driver_genes = set()

with open('data/COSMIC_Census_allMon May  7 17_11_47 2018.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')
    reader.next() # skip header
    for row in reader:
        driver_genes.add(row[0])



vogelstein_drivers = set(
['ABL1','ACVR1B','AKT1','ALK','APC','AR','ARID1A','ARID1B',
'ARID2','ASXL1','ATM','ATRX','AXIN1','B2M','BAP1','BCL2',
'BCOR','BRAF','BRCA1','BRCA2','CARD11','CASP8','CBL','CDC73',
'CDH1','CDKN2A','CEBPA','CIC','CREBBP','CRLF2','CSF1R','CTNNB1',
'CYLD','DAXX','DNMT1','DNMT3A','EGFR','EP300','ERBB2','EZH2',
'FAM123B','FBXW7','FGFR2','FGFR3','FLT3','FOXL2','FUBP1','GATA1','GATA2','GATA3',
'GNA11','GNAQ','GNAS','H3F3A','HIST1H3B','HNF1A','HRAS','IDH1','IDH2',
'JAK1','JAK2','JAK3','KDM5C','KDM6A','KIT','KLF4','KRAS','MAP2K1','MAP3K1',
'MED12','MEN1','MET','MLH1','MLL2','MLL3','MPL','MSH2',
'MSH6','MYD88','NCOR1','NF1','NF2','NFE2L2','NOTCH1','NOTCH2',
'NPM1','NRAS','PAX5','PBRM1','PDGFRA','PHF6','PIK3CA','PIK3R1',
'PPP2R1A','PRDM1','PTCH1','PTEN','PTPN11','RB1','RET','RNF43',
'RUNX1','SETD2','SETBP1','SF3B1','SMAD2','SMAD4','SMARCA4','SMARCB1',
'SMO','SOCS1','SOX9','SPOP','SRSF2','STAG2','STK11','TET2',
'TNFAIP3','TRAF7','TP53','TSC1','TSHR','U2AF1','VHL','WT1']

)

# Check how many of the driver genes contain frequent kmers

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

print "Found %d driver genes from COSMIC database (https://cancer.sanger.ac.uk/census)" % \
    len(driver_genes)

driver_that_contain_frequent_kmers = set()
for rec in SeqIO.parse(open(human_proteome), 'fasta'):
    seq = str(rec.seq)
    uniqueIdentifier, entryName, proteinName, organismName, geneName = \
        utils.parse_UniProtKB_header(rec.description)
    if geneName in driver_genes and \
        any([kmer in seq for kmer in most_frequent_kmers]):
        driver_that_contain_frequent_kmers.add(geneName)

print "%d of them contain one of the top-50 frequent kmers" % len(driver_that_contain_frequent_kmers)

print

print "Found %d driver genes from Vogelstein database" % \
    len(vogelstein_drivers)

driver_that_contain_frequent_kmers = set()
for rec in SeqIO.parse(open(human_proteome), 'fasta'):
    seq = str(rec.seq)
    uniqueIdentifier, entryName, proteinName, organismName, geneName = \
        utils.parse_UniProtKB_header(rec.description)
    if geneName in vogelstein_drivers and \
        any([kmer in seq for kmer in most_frequent_kmers]):
        driver_that_contain_frequent_kmers.add(geneName)

print "%d of them contain one of the top-50 frequent kmers" % len(driver_that_contain_frequent_kmers)

print


genes_that_contain_frequent_kmers = set()
for rec in SeqIO.parse(open(human_proteome), 'fasta'):
    seq = str(rec.seq)
    uniqueIdentifier, entryName, proteinName, organismName, geneName = \
        utils.parse_UniProtKB_header(rec.description)
    if any([kmer in seq for kmer in most_frequent_kmers]):
        genes_that_contain_frequent_kmers.add(geneName)

print "In total, %d genes contain one of the top-50 kmers" % len(genes_that_contain_frequent_kmers)

print "==="

for freq_kmer in most_frequent_kmers:
    proteins_that_contain_fequent_kmer = set()
    driver_that_contain_frequent_kmer = set()
    for rec in SeqIO.parse(open(human_proteome), 'fasta'):
        seq = str(rec.seq)
        uniqueIdentifier, entryName, proteinName, organismName, geneName = \
            utils.parse_UniProtKB_header(rec.description)
        if geneName in driver_genes and freq_kmer in seq:
            driver_that_contain_frequent_kmer.add(geneName)
        if freq_kmer in seq:
            proteins_that_contain_fequent_kmer.add(geneName)
    print "%d proteins contain the kmer %s, %d of them are driver genes." % \
          (len(proteins_that_contain_fequent_kmer),
           freq_kmer,
           len(driver_that_contain_frequent_kmer))