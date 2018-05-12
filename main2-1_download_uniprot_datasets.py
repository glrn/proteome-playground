import csv
import urllib

EYKARYOTES_ID_FILE  = "Organism_UniProt_IDs_Eukaryotes.csv"
BACTERIA_ID_FILE    = "Organism_UniProt_IDs_Bacteria.csv"
UNIPROT_ORGANISM_FASTA_LINK = 'http://www.uniprot.org/uniprot/?query=organism:%d&format=fasta'

with open(EYKARYOTES_ID_FILE, 'rb') as csvfile:
    records = csv.DictReader(csvfile)
    for row in records:
        organism_name = row['Scientific']
        organism_id = int(row['UniProt Organism ID'])
        uniprot_fasta_url = UNIPROT_ORGANISM_FASTA_LINK % organism_id
        dest_file = "data/Uniprot Proteomes/Eukaryotes/%s.fasta" % organism_name
        print 'Downloading FASTA file of %s from UniProt --> %s' % (organism_name, dest_file)
        retriever = urllib.URLopener()
        retriever.retrieve(uniprot_fasta_url, dest_file)

with open(BACTERIA_ID_FILE, 'rb') as csvfile:
    records = csv.DictReader(csvfile)
    for row in records:
        organism_name = row['Scientific']
        organism_id = int(row['UniProt Organism ID'])
        uniprot_fasta_url = UNIPROT_ORGANISM_FASTA_LINK % organism_id
        dest_file = "data/Uniprot Proteomes/Bacteria/%s.fasta" % organism_name
        print 'Downloading FASTA file of %s from UniProt --> %s' % (organism_name, dest_file)
        retriever = urllib.URLopener()
        retriever.retrieve(uniprot_fasta_url, dest_file)
