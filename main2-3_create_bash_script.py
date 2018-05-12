from os import listdir
from os.path import isfile, join, splitext

data_path = "data/Uniprot Proteomes/Eukaryotes"
fastafiles = [data_path + '/' + f for f in listdir(data_path) if isfile(join(data_path, f)) and
             splitext(join(data_path, f))[1] == '.fasta']

for file in fastafiles:
    #print 'python main2-2_kmer_counting.py "{}" > /dev/null &'.format(file)
    print 'python main2-2_kmer_counting.py "{}" --similar-diluted > /dev/null &'.format(file)