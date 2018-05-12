from Bio import SeqIO
import sys
import urllib2
import json
# local imports
from progressbar import Progressbar
import params
import utils
import random

THE_TWENTY_AMINO_ACIDS = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

def random_human_genes(num):
    genes_set = set()
    for rec in SeqIO.parse(open("data\UniProt full Proteomes\uniprot-proteome-full Homo sapiens (71599 proteins).fasta"), 'fasta'):
        seq = str(rec.seq)
        uniqueIdentifier, entryName, proteinName, organismName, geneName = \
            utils.parse_UniProtKB_header(rec.description)
        genes_set.add(geneName)
    return random.sample(genes_set, num)

def get_human_proteins_that_contain_saar(SAAR):
    rcp_set = set()
    genes_set = set()
    for rec in SeqIO.parse(open("data\UniProt full Proteomes\uniprot-proteome-full Homo sapiens (71599 proteins).fasta"), 'fasta'):
        seq = str(rec.seq)
        uniqueIdentifier, entryName, proteinName, organismName, geneName = \
            utils.parse_UniProtKB_header(rec.description)
        genes_set.add(geneName)
        if SAAR in seq:
            rcp_set.add(geneName)
    return rcp_set

def get_enrichment(genes_set):
    string_api_url = "https://string-db.org/api"
    output_format = "json"
    method = "enrichment"

    my_genes = genes_set

    species = "9606"
    my_app = "www.awesome_app.org"

    ## Construct the request

    request_url = string_api_url + "/" + output_format + "/" + method + "?"
    request_url += "identifiers=" + "%0d".join(my_genes)
    request_url += "&" + "species=" + species
    request_url += "&" + "caller_identity=" + my_app

    ## Call STRING

    try:
        response = urllib2.urlopen(request_url)
    except urllib2.HTTPError as err:
        error_message = err.read()
        print error_message
        sys.exit()

    ## Read and parse the results

    result = response.readline()

    if result:
        data = json.loads(result)
        for row in data:
            term = row["term"]
            preferred_names = ",".join(row["preferredNames"])
            fdr = row["fdr"]
            description = row["description"]

            if fdr < 0.01:
                print "\t".join([term, preferred_names, str(fdr), description])

def get_ppi_enrichment(genes_set):
    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "ppi_enrichment"

    my_genes = genes_set

    species = "9606"
    my_app = "www.awesome_app.org"

    ## Construct the request

    request_url = string_api_url + "/" + output_format + "/" + method + "?"
    request_url += "identifiers=" + "%0d".join(my_genes)
    request_url += "&" + "species=" + species
    request_url += "&" + "caller_identity=" + my_app

    ## Call STRING

    try:
        response = urllib2.urlopen(request_url)
    except urllib2.HTTPError as err:
        error_message = err.read()
        print error_message
        sys.exit()

    ## Read and parse the results

    result = response.readline()

    if result:
        pvalue = result.strip().split("\t")[5]

    return pvalue

"""
genes = ['ZFR', 'NIPA1', 'KMT5A', 'NKX2-3', 'NKX2-4', 'OTUD7A', 'LHX2', 'BHLHE41', 'MAP3K4', 'ZSWIM6', 'SLC12A2', 'PTBP2', 'IRX1', 'SP8', 'SP5', 'IRF2BPL', 'EVX2', 'RPL14', 'GATA6', 'NKX6-1', 'RBM47', 'ARID1B', 'MAZ', 'RUNX2', 'PABPN1', 'MNX1', 'CHD3', 'EGR2', 'ASXL2', 'ZBTB8B', 'EN1', 'ZIC3', 'NLK', 'GDF11', 'PRDM12', 'HOXD11', 'FBRS', 'RBM24', 'CASKIN1', 'ALKBH5', 'ZNF358', 'LCORL', 'BTBD2', 'KDM1A', 'TLE4', 'TBX2', 'TMEM200C', 'LHFPL3', 'SKIDA1', 'HOXA13', 'PRPF40A', 'REPS2', 'FOXB2', 'CLASRP', 'KLHL11', 'FBXO41', 'ASCL1', 'FOXD3', 'FEV', 'FOXL2', 'PHLPP1', 'OLIG2', 'TBL1XR1', 'IRS4', 'ID4', 'PRKD1', 'SKOR1', 'HMX3', 'RBM4', 'MSH3', 'JUND', 'HAND2', 'ELN', 'TMEM158', 'CDK13', 'EOMES', 'POU3F1', 'POU4F1', 'PHOX2B', 'ZNF804A', 'SOX21', 'ARX', 'PRDM8', 'DMRTA1', 'CCDC177', 'GSX2']
genes = ['TIMELESS', 'IL27', 'NCBP3', 'PRKCSH', 'EHMT2', 'TSPOAP1', 'BCL11A', 'ZFHX3', 'HNRNPU', 'GTF3C3', 'ANP32E', 'EP400', 'CACNA1F', 'SLC4A3', 'AEBP2', 'SLC24A1', 'ARID4B', 'CSRNP3', 'FKBP8', 'C11orf95', 'DAXX', 'CHIC1', 'CERS3', 'ZBTB47', 'CNGB1', 'EIF5B', 'EIF3D', 'ALMS1', 'ADRA2B', 'CNKSR2', 'HMGB3', 'CNPY4', 'CUL9', 'WDR70', 'KDM2A', 'KCNA4', 'CLSTN1', 'PODN', 'TTBK1', 'MYT1', 'SHROOM4', 'WEE1', 'PODXL2', 'NGEF', 'STAC3', 'PLPPR3', 'RTN4', 'TRIM26', 'SPRYD3', 'ACIN1', 'SENP3', 'CELSR2', 'CDK11B', 'CENPB', 'PELP1', 'SMYD5', 'BAZ1B', 'PPARGC1B', 'SCAF1', 'SNAPC5', 'MYO15B', 'ATRX', 'PAXIP1', 'HRC', 'RBM19', 'PHF23', 'NEUROD2', 'YTHDC1', 'CCDC136', 'TAOK2', 'PLEKHG5', 'SMARCA4', 'FAM9A', 'WIZ', 'SETD1B', 'VSIG10', 'FBXO3', 'NCL', 'SKIDA1', 'RRP12', 'KAT6B', 'DCAF8', 'TPRN', 'SLC4A1AP', 'RANGAP1']
genes = ['RANBP9', 'HS3ST4', 'BAG6', 'ESPN', 'FNBP4', 'ZNF318', 'FBXO11', 'ZFHX3', 'SF3B2', 'SETD1A', 'SMN1', 'C11orf95', 'PEG10', 'SEMA4G', 'PXK', 'JMY', 'MYPOP', 'WASF2', 'SIK2', 'ASXL3', 'FMNL1', 'ZNF827', 'MAP3K4', 'AHDC1', 'RAPH1', 'ZNF341', 'ZIC5', 'PCLO', 'KDM6B', 'VIRMA', 'PRIMA1', 'FOXB2', 'DOT1L', 'GJD3', 'ADGRB1', 'TRERF1', 'MICALCL', 'DIAPH3', 'RNF44', 'DAB2IP', 'R3HDM1', 'HOXB4', 'SP5', 'ZNF462', 'TRRAP', 'SALL2', 'ALG13', 'WAS', 'MEX3C', 'SRPK2', 'FHDC1', 'IRS4', 'ZC3H4', 'USP51', 'PRR12', 'WASH2P', 'HTT', 'WIPF3', 'ANTXRL', 'PPP3CB', 'SKOR2', 'LMOD2']
genes = ['CIR1', 'PPRC1', 'PHLPP2', 'MPRIP', 'PRG4', 'EDC4', 'RANBP10', 'NR1D1', 'ATXN7L1', 'PRDM2', 'BHMG1', 'JMJD6', 'CHD8', 'MLLT6', 'ZBTB4', 'KCNMA1', 'MLLT3', 'POU6F2', 'ZNF865', 'FAM193B', 'LZTS2', 'JCAD', 'EPHB6', 'HEG1', 'TPRXL', 'SALL1', 'TSHZ3', 'NOLC1', 'ARL6IP4', 'GTF2IRD1', 'ATN1', 'SHE', 'SALL2', 'DACH1', 'SERTAD4', 'MEX3B', 'BRD2', 'NKX6-1', 'NAF1', 'SETD1A', 'E2F4', 'MBTPS2', 'HDGFL2', 'RLIM', 'SRRM2', 'TMEM40', 'PRR13', 'DCAF5', 'PNN', 'PPIG', 'RNF126', 'BCL6B', 'SYN2', 'NFE2L1', 'TNRC18']
genes = ['HMGB3', 'NISCH', 'PPP4R2', 'BCL11B', 'DIEXF', 'KCNA4', 'APLP2', 'TTBK1', 'ZBTB4', 'PPARGC1B', 'GOLGA2', 'VGLL3', 'MYT1', 'PLPPR3', 'BAHCC1', 'RANGAP1', 'PELP1', 'HNRNPUL2', 'PTMS', 'KAT6B', 'PDE4C', 'ZNF428', 'PAK3']
genes = ['ZNF215', 'ZNF536', 'ZNF554', 'ZNF319', 'ZSCAN29', 'ZNF233', 'ZNF449', 'ZNF597', 'ZNF75D', 'ZNF781', 'ZNF18', 'ZNF768', 'ZNF236', 'KLF8', 'ZNF394', 'ZNF429']
print get_ppi_enrichment(genes)
sys.exit()
"""

print 'PPI enrichment value (according to String): %s' % get_enrichment(['FOXP2', 'LIN7A', 'ZFHX3', 'SATB1', 'SIK3', 'CELF3', 'TMEM229A', 'POLG', 'GLG1', 'EP400', 'DCP1B', 'DENND4B', 'SMARCA2', 'THAP11', 'PHLDA1', 'NCOR2', 'NCOA3', 'ZNF384', 'MAML2', 'FAM157A', 'ABCF1', 'USF3', 'CACNA1A', 'TFEB', 'MEF2A', 'POU6F2', 'RAI1', 'MED12', 'PHC1', 'MED15', 'MAGI1', 'NUMBL', 'POU3F2', 'KMT2D', 'TSC22D1', 'MN1', 'ASCL1', 'R3HDM2', 'IRF2BPL', 'KCNN3', 'TFCP2', 'AAK1', 'AR', 'ATN1', 'CHERP', 'FRMPD3', 'NCOA6', 'ST6GALNAC5', 'MAML3', 'NFAT5', 'ARID3B', 'FAM171B', 'VEZF1', 'ARID1B', 'RUNX2', 'MAMLD1', 'HTT', 'TOX3', 'ATXN8', 'TBP', 'ATXN7', 'ATXN2', 'ATXN3', 'CREBBP', 'ATXN1'])
sys.exit()

for AA in THE_TWENTY_AMINO_ACIDS:
    rcp_set = get_human_proteins_that_contain_saar(AA*10)
    print 'Total number of human genes containing %s: %d' % (AA*10, len(rcp_set))
    if len(rcp_set) > 0:
        print 'PPI enrichment value (according to String): %s' % get_ppi_enrichment(rcp_set)
        print 'PPI enrichment value (according to String): %s' % get_enrichment(rcp_set)

    print '-----'
