
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
