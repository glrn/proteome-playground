
def parse_UniProtKB_header(fastaHeader):
    # >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName[ GN=GeneName]PE=ProteinExistence SV=SequenceVersion
    fastaHeader = fastaHeader[fastaHeader.find('|')+1:]
    uniqueIdentifier = fastaHeader[:fastaHeader.find('|')]
    fastaHeader = fastaHeader[fastaHeader.find('|') + 1:]
    entryName = fastaHeader[:fastaHeader.find(' ')]
    fastaHeader = fastaHeader[fastaHeader.find(' ') + 1:]
    proteinName = fastaHeader[:fastaHeader.find('OS=')-1]
    fastaHeader = fastaHeader[fastaHeader.find('OS='):]
    organismName = fastaHeader[fastaHeader.find('OS=')+3:fastaHeader.find('GN=')-1]
    fastaHeader = fastaHeader[fastaHeader.find(organismName)+len(organismName)+1:]
    geneName = fastaHeader[fastaHeader.find('GN=')+3:fastaHeader.find(' ')]
    return uniqueIdentifier, entryName, proteinName, organismName, geneName

def parse_NCBIRefSeq_header(fastaHeader):
    organismName = fastaHeader[fastaHeader.rfind('[')+1:fastaHeader.rfind(']')]
    geneName = fastaHeader[fastaHeader.find(' ')+1:fastaHeader.rfind('[')]
    return organismName, geneName

dissimilarity_mapping = {} # dissimilarity_mapping((protNameA,protNameB)) = TRUE/FALSE

def proteins_are_dissimilar(protNameA, protNameB, protA, protB):
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    if (protNameA, protNameB) in dissimilarity_mapping.keys():
        return dissimilarity_mapping[(protNameA, protNameB)]
    score = pairwise2.align.globalxs(protA, protB, -0.5, -0.1, score_only=True)
    smallerProt = min(len(protA),len(protB))
    if score < smallerProt * 0.4:
        dissimilarity_mapping[(protNameA, protNameB)] = True
        return True # proteins are dissimilar enough
    else:
        dissimilarity_mapping[(protNameA, protNameB)] = False
        return False


def hamdist(str1, str2):
    """Count the # of differences between equal length strings str1 and str2"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs