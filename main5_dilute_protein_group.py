from Bio import SeqIO
# local imports
from progressbar import Progressbar
import params
import utils

class Protein:
    """
    Hold information of gene name, amino acids, and all k-mers of amino acids.
    """

    def __init__(self, geneName, sequence, k):
        self.geneName = geneName
        self.seq = sequence


all_proteins = list() # all_proteins[i] = PROTEIN()_OBJECT
protein_seq = dict() # protein_seq[GENE_NAME] = AMINO_ACID_SEQUENCE

print "Reading Uniprot file..."
created_protein_names = set() # prevent creation of two similar protein objects
for rec in SeqIO.parse(open(params.PROTEOME_FILE), 'fasta'):
    seq = str(rec.seq)
    uniqueIdentifier, entryName, proteinName, organismName, geneName = \
        utils.parse_UniProtKB_header(rec.description)
    if geneName == '':
        print 'Ignoring unknown gene: %s' % rec.description
        continue
    if geneName in created_protein_names:
        print 'Ignoring duplicate gene: %s' % rec.description
        continue
    # create a new Protein object
    created_protein_names.add(geneName)
    protein_seq[geneName] = seq
    all_proteins.append(Protein(geneName, seq, params.K))

g1 = ['FOXP2', 'EP400', 'SMARCA2', 'THAP11', 'NCOA3', 'MAML2', 'MAML3', 'BMP2K',
     'MED12', 'MED15', 'POU3F2', 'KMT2D', 'MN1', 'IRF2BPL', 'AR', 'KCNN3',
     'FRMPD3', 'NCOA6', 'NFAT5', 'TBP', 'RUNX2', 'MAMLD1', 'HTT', 'ATXN8',
     'ATXN2', 'ATXN1']
g2 = ['ARL6IP4', 'MBTPS2', 'TPRXL', 'TNRC18', 'PPRC1', 'KCNMA1', 'MLLT3',
     'SETD1A', 'ZBTB4', 'ZNF865', 'DACH1', 'TMEM40', 'SRRM2']
g3 = ['ZIC5', 'FNBP4', 'ZFHX4', 'FMNL1', 'SETD1A', 'RAPH1', 'HTT', 'PCLO',
     'ALG13', 'SKOR2', 'LMOD2']
g4 = ['EHMT2', 'ZBTB47', 'CDK11B', 'PELP1', 'CHIC1', 'SCAF1', 'DCAF8L2', 'KAT6B',
     'TTBK1', 'MYT1']
g5 = ['POU3F3', 'SOX21', 'ARX', 'SP8', 'MAZ', 'HOXA13', 'IRF2BPL', 'EVX2',
     'RBM24', 'FOXB2']
g6 = ['POU3F2', 'CAPNS1', 'FUS', 'HOXB3', 'NPAS3', 'AR', 'DACH1', 'ZNF503']
g7 = ['ZNF334', 'ZNF619', 'ZKSCAN5', 'GLI4', 'ZNF473', 'ZSCAN20']
g8 = ['ZNF263', 'ZNF487', 'ZNF268', 'PRDM9', 'ZNF629']
g9 = ['ZNF840P', 'ZSCAN31', 'ZFP62', 'GLI4', 'ZNF561']
g10 = ['ZNF487', 'ZNF629', 'ZNF202', 'PRDM9', 'ZNF577']

all_prots = g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + g9 + g10
genes = set()
genes.update(set(['ZFR', 'NIPA1', 'KMT5A', 'NKX2-3', 'NKX2-4', 'OTUD7A', 'LHX2', 'BHLHE41', 'MAP3K4', 'ZSWIM6', 'SLC12A2', 'PTBP2', 'IRX1', 'SP8', 'SP5', 'IRF2BPL', 'EVX2', 'RPL14', 'GATA6', 'NKX6-1', 'RBM47', 'ARID1B', 'MAZ', 'RUNX2', 'PABPN1', 'MNX1', 'CHD3', 'EGR2', 'ASXL2', 'ZBTB8B', 'EN1', 'ZIC3', 'NLK', 'GDF11', 'PRDM12', 'HOXD11', 'FBRS', 'RBM24', 'CASKIN1', 'ALKBH5', 'ZNF358', 'LCORL', 'BTBD2', 'KDM1A', 'TLE4', 'TBX2', 'TMEM200C', 'LHFPL3', 'SKIDA1', 'HOXA13', 'PRPF40A', 'REPS2', 'FOXB2', 'CLASRP', 'KLHL11', 'FBXO41', 'ASCL1', 'FOXD3', 'FEV', 'FOXL2', 'PHLPP1', 'OLIG2', 'TBL1XR1', 'IRS4', 'ID4', 'PRKD1', 'SKOR1', 'HMX3', 'RBM4', 'MSH3', 'JUND', 'HAND2', 'ELN', 'TMEM158', 'CDK13', 'EOMES', 'POU3F1', 'POU4F1', 'PHOX2B', 'ZNF804A', 'SOX21', 'ARX', 'PRDM8', 'DMRTA1', 'CCDC177', 'GSX2']))
genes.update(set(['TIMELESS', 'IL27', 'NCBP3', 'PRKCSH', 'EHMT2', 'TSPOAP1', 'BCL11A', 'ZFHX3', 'HNRNPU', 'GTF3C3', 'ANP32E', 'EP400', 'CACNA1F', 'SLC4A3', 'AEBP2', 'SLC24A1', 'ARID4B', 'CSRNP3', 'FKBP8', 'C11orf95', 'DAXX', 'CHIC1', 'CERS3', 'ZBTB47', 'CNGB1', 'EIF5B', 'EIF3D', 'ALMS1', 'ADRA2B', 'CNKSR2', 'HMGB3', 'CNPY4', 'CUL9', 'WDR70', 'KDM2A', 'KCNA4', 'CLSTN1', 'PODN', 'TTBK1', 'MYT1', 'SHROOM4', 'WEE1', 'PODXL2', 'NGEF', 'STAC3', 'PLPPR3', 'RTN4', 'TRIM26', 'SPRYD3', 'ACIN1', 'SENP3', 'CELSR2', 'CDK11B', 'CENPB', 'PELP1', 'SMYD5', 'BAZ1B', 'PPARGC1B', 'SCAF1', 'SNAPC5', 'MYO15B', 'ATRX', 'PAXIP1', 'HRC', 'RBM19', 'PHF23', 'NEUROD2', 'YTHDC1', 'CCDC136', 'TAOK2', 'PLEKHG5', 'SMARCA4', 'FAM9A', 'WIZ', 'SETD1B', 'VSIG10', 'FBXO3', 'NCL', 'SKIDA1', 'RRP12', 'KAT6B', 'DCAF8', 'TPRN', 'SLC4A1AP', 'RANGAP1']))
genes.update(set(['FOXP2', 'LIN7A', 'ZFHX3', 'SATB1', 'SIK3', 'CELF3', 'TMEM229A', 'POLG', 'GLG1', 'EP400', 'DCP1B', 'DENND4B', 'SMARCA2', 'THAP11', 'PHLDA1', 'NCOR2', 'NCOA3', 'ZNF384', 'MAML2', 'FAM157A', 'ABCF1', 'USF3', 'CACNA1A', 'TFEB', 'MEF2A', 'POU6F2', 'RAI1', 'MED12', 'PHC1', 'MED15', 'MAGI1', 'NUMBL', 'POU3F2', 'KMT2D', 'TSC22D1', 'MN1', 'ASCL1', 'R3HDM2', 'IRF2BPL', 'KCNN3', 'TFCP2', 'AAK1', 'AR', 'ATN1', 'CHERP', 'FRMPD3', 'NCOA6', 'ST6GALNAC5', 'MAML3', 'NFAT5', 'ARID3B', 'FAM171B', 'VEZF1', 'ARID1B', 'RUNX2', 'MAMLD1', 'HTT', 'TOX3', 'ATXN8', 'TBP', 'ATXN7', 'ATXN2', 'ATXN3', 'CREBBP', 'ATXN1']))
genes.update(set(['RANBP9', 'HS3ST4', 'BAG6', 'ESPN', 'FNBP4', 'ZNF318', 'FBXO11', 'ZFHX3', 'SF3B2', 'SETD1A', 'SMN1', 'C11orf95', 'PEG10', 'SEMA4G', 'PXK', 'JMY', 'MYPOP', 'WASF2', 'SIK2', 'ASXL3', 'FMNL1', 'ZNF827', 'MAP3K4', 'AHDC1', 'RAPH1', 'ZNF341', 'ZIC5', 'PCLO', 'KDM6B', 'VIRMA', 'PRIMA1', 'FOXB2', 'DOT1L', 'GJD3', 'ADGRB1', 'TRERF1', 'MICALCL', 'DIAPH3', 'RNF44', 'DAB2IP', 'R3HDM1', 'HOXB4', 'SP5', 'ZNF462', 'TRRAP', 'SALL2', 'ALG13', 'WAS', 'MEX3C', 'SRPK2', 'FHDC1', 'IRS4', 'ZC3H4', 'USP51', 'PRR12', 'WASH2P', 'HTT', 'WIPF3', 'ANTXRL', 'PPP3CB', 'SKOR2', 'LMOD2']))
genes.update(set(['CIR1', 'PPRC1', 'PHLPP2', 'MPRIP', 'PRG4', 'EDC4', 'RANBP10', 'NR1D1', 'ATXN7L1', 'PRDM2', 'BHMG1', 'JMJD6', 'CHD8', 'MLLT6', 'ZBTB4', 'KCNMA1', 'MLLT3', 'POU6F2', 'ZNF865', 'FAM193B', 'LZTS2', 'JCAD', 'EPHB6', 'HEG1', 'TPRXL', 'SALL1', 'TSHZ3', 'NOLC1', 'ARL6IP4', 'GTF2IRD1', 'ATN1', 'SHE', 'SALL2', 'DACH1', 'SERTAD4', 'MEX3B', 'BRD2', 'NKX6-1', 'NAF1', 'SETD1A', 'E2F4', 'MBTPS2', 'HDGFL2', 'RLIM', 'SRRM2', 'TMEM40', 'PRR13', 'DCAF5', 'PNN', 'PPIG', 'RNF126', 'BCL6B', 'SYN2', 'NFE2L1', 'TNRC18']))
genes.update(set(['CAPNS1', 'FUS', 'PI4KA', 'NPAS3', 'FUBP1', 'ZFHX3', 'RFX1', 'LURAP1L', 'SRSF11', 'AEBP2', 'ANKRD17', 'XYLT1', 'SOX1', 'DHX36', 'FEZF2', 'U2AF1', 'GDF7', 'WDR26', 'FOXF1', 'FBXL17', 'ONECUT2', 'TRA2A', 'SP8', 'POU3F2', 'TAL1', 'HOXD9', 'ARIH1', 'HOXB3', 'PITPNM2', 'EMC10', 'HOXC13', 'AR', 'CCDC6', 'DACH1', 'MEX3D', 'EVX2', 'ZNF503', 'UNC80', 'HCN1', 'ARID1B', 'GSK3A', 'POU4F2', 'TCF7L1', 'YEATS2', 'PURA', 'SHOX2', 'MAF', 'WIPF1', 'FOXC1', 'FBL']))
genes.update(set(['GTF3C5', 'HMGB3', 'BCL11A', 'PELP1', 'GTF3C3', 'ZNF830', 'SLC4A3', 'KCNA4', 'DAXX', 'CHIC1', 'CNGB1', 'SUPT5H', 'EIF3D', 'HNRNPU', 'PPARGC1B', 'CSRNP3', 'RTN4', 'PIK3R5', 'KIF21A', 'MYT1', 'ZBTB7C', 'PLPPR3', 'HOMEZ', 'MIER1', 'TRIM26', 'RANGAP1', 'BAZ1B', 'TTBK1', 'ZEB1', 'MAPK8IP2', 'NCL', 'PIAS4', 'ZBTB47', 'KAT6B']))
genes.update(set(['NCBP3', 'HMGB3', 'ZFHX3', 'PELP1', 'FBXO3', 'GTF3C5', 'ARID4B', 'EIF5B', 'EHMT2', 'VSIG10', 'FKBP8', 'NAP1L2', 'CSRNP3', 'RTN4', 'GZF1', 'KDM2A', 'KIF21A', 'MYT1', 'PRDM2', 'ZBTB7C', 'PLPPR3', 'HOMEZ', 'ACIN1', 'SENP3', 'BAZ1B', 'SCAF1', 'ATRX', 'NEUROD1', 'PCGF6', 'ANP32E', 'SPRYD3', 'PIAS4', 'KAT6B', 'WIZ']))
genes.update(set(['HMX3', 'TLE3', 'BHLHE22', 'SP5', 'TBX2', 'NKX2-3', 'PITX3', 'EGR2', 'HOXA13', 'PBX3', 'EOMES', 'ZSWIM6', 'SPTLC2', 'FBXL17', 'ZNF703', 'ZIC2', 'NLK', 'SP8', 'FOXD3', 'IRF2BPL', 'EVX2', 'GATA6', 'FOXB2', 'NKX6-1', 'SOX21', 'RBM47', 'ID4', 'MAZ', 'GSX2', 'RBM24', 'PRDM8', 'ZNF395']))
genes.update(set(['NIPA1', 'JUND', 'TBX2', 'ZNF358', 'MNX1', 'SKIDA1', 'HOXA13', 'REPS2', 'ZSWIM6', 'ZIC5', 'SPTLC2', 'SLC12A2', 'PTBP2', 'FOXB2', 'GDF11', 'PABPN1', 'HOXD8', 'POU4F1', 'FOXL2', 'C6orf223', 'PHOX2B', 'ZNF804A', 'NRG3', 'SOX21', 'DMRTA2', 'ARID1B', 'SOWAHA', 'PRDM8', 'ZNF395', 'SKOR1']))
genes.update(set(['HS3ST4', 'ESPN', 'ONECUT3', 'ENAH', 'FBXO11', 'FMN1', 'ZFHX3', 'INF2', 'JMY', 'CDK12', 'ASXL3', 'RAPH1', 'ZNF341', 'PCLO', 'KDM6B', 'DOT1L', 'TSPYL2', 'DIAPH3', 'R3HDM1', 'FASLG', 'UBR4', 'MEX3C', 'SRPK2', 'SHTN1', 'FMNL2', 'ANTXRL', 'SKOR2', 'LMOD2']))
genes.update(set(['FOXK1', 'WASF1', 'ENAH', 'FBXO11', 'FMN1', 'ZFHX3', 'INF2', 'JMY', 'ESPN', 'FMNL2', 'ASXL3', 'YLPM1', 'RAPH1', 'ZNF341', 'PCLO', 'KDM6B', 'DOT1L', 'DIAPH3', 'R3HDM1', 'FASLG', 'MEX3C', 'SRPK2', 'NAF1', 'FHDC1', 'SHTN1', 'PRR11', 'ANTXRL', 'LMOD2']))
genes.update(set(['ZNF316', 'HMGB2', 'NISCH', 'ZFHX3', 'FBXO3', 'PIAS4', 'TSPOAP1', 'ARID4B', 'EIF5B', 'EHMT2', 'CACNA1D', 'CLIP3', 'GZF1', 'MYT1', 'KRI1', 'PLPPR3', 'RTN4', 'SENP3', 'SCAF1', 'ZNF777', 'ATRX', 'TTBK1', 'BOD1L1', 'ANP32E', 'PRDM2', 'KAT6B']))
genes.update(set(['CPEB3', 'HAND2', 'PITX3', 'SOX1', 'ZNF703', 'CDK13', 'ASXL2', 'HOXA13', 'TMEM200C', 'FBXL17', 'FOXB2', 'IRX1', 'CLASRP', 'SP8', 'KLHL11', 'FBXO41', 'EVX2', 'GATA6', 'PBX3', 'NKX6-1', 'TBL1XR1', 'IRS4', 'ARID1A', 'CCDC177', 'TOX4']))
genes.update(set(['ZNF830', 'SLC4A3', 'BCORL1', 'C11orf95', 'KCNA4', 'SKIDA1', 'HNRNPU', 'CACNA1F', 'CNKSR2', 'CNPY4', 'MAPK8IP2', 'WEE1', 'PRX', 'TPRN', 'MIER1', 'CENPB', 'CDK11B', 'SMYD5', 'TAOK2', 'NEUROD2', 'PLEKHG5', 'PRM3', 'ZBTB47', 'KAT6B', 'ZNF428']))
genes.update(set(['TRERF1', 'ZNF318', 'FMN1', 'PELP1', 'ACR', 'DMRTB1', 'SETD1B', 'C11orf95', 'MYPOP', 'ASXL3', 'FMNL1', 'PXK', 'RAPH1', 'ZIC5', 'PRIMA1', 'WHAMM', 'KDM6B', 'ESPN', 'GJD3', 'SCAF11', 'RNF44', 'SP5', 'ALG13', 'LMOD2']))
genes.update(set(['DTX3', 'FMN1', 'SETBP1', 'ZBTB12', 'WASF1', 'ASXL3', 'YLPM1', 'CAP2', 'ZNF341', 'PCLO', 'VCPIP1', 'GRID2IP', 'DOT1L', 'ARHGEF11', 'FASLG', 'MTFR1', 'MEX3C', 'NAF1', 'SHTN1', 'PRR11', 'WIPF3', 'ANTXRL', 'PDE4D']))
genes.update(set(['HMGB3', 'NISCH', 'PPP4R2', 'BCL11B', 'DIEXF', 'KCNA4', 'APLP2', 'TTBK1', 'ZBTB4', 'PPARGC1B', 'GOLGA2', 'VGLL3', 'MYT1', 'PLPPR3', 'BAHCC1', 'RANGAP1', 'PELP1', 'HNRNPUL2', 'PTMS', 'KAT6B', 'PDE4C', 'ZNF428', 'PAK3']))
genes.update(set(['FOXK1', 'ZFHX4', 'FMN1', 'SETBP1', 'WASF1', 'ASXL3', 'YLPM1', 'ZNF341', 'PCLO', 'VCPIP1', 'GRID2IP', 'DOT1L', 'DAAM2', 'R3HDM1', 'FASLG', 'MTFR1', 'MEX3C', 'NAF1', 'SHTN1', 'PRR11', 'WIPF3', 'ANTXRL']))
genes.update(set(['HMGB3', 'BCL11B', 'ZNF830', 'GTF3C5', 'KCNA4', 'CHIC1', 'CNGB1', 'SUPT5H', 'EIF3D', 'PPARGC1B', 'TMED8', 'KIF21A', 'MYT1', 'ZBTB7C', 'PLPPR3', 'RTN4', 'MIER1', 'RANGAP1', 'TTBK1', 'MAPK8IP2', 'PIAS4', 'KAT6B']))
genes.update(set(['ANO8', 'PPP6R1', 'PELP1', 'DIEXF', 'KCNA4', 'APLP2', 'ZBTB4', 'GOLGA2', 'RBM10', 'MYT1', 'PLPPR3', 'BAHCC1', 'RANGAP1', 'SCAF1', 'SALL2', 'TTBK1', 'FAM193A', 'TAF1', 'PRDM2', 'KAT6B', 'ITGAE', 'ZNF428']))
genes.update(set(['HMGB3', 'NISCH', 'BCL11B', 'DIEXF', 'GTF3C5', 'APLP2', 'CWC22', 'KCNA4', 'CHIC1', 'PPARGC1B', 'GOLGA2', 'VGLL3', 'MYT1', 'ZBTB7C', 'PLPPR3', 'MIER1', 'RANGAP1', 'TTBK1', 'MBD3', 'KAT6B', 'ZNF428', 'PAK3']))
genes.update(set(['CEBPB', 'PXK', 'BAG6', 'ZIC5', 'FNBP4', 'FIGN', 'PRR12', 'WASH2P', 'DAB2IP', 'TRRAP', 'USP51', 'AHDC1', 'IRS4', 'SETD1B', 'RANBP9', 'SCX', 'GABBR2', 'ALG13', 'GJD3', 'HNRNPUL1', 'SEMA4G']))
genes.update(set(['EHMT2', 'TULP1', 'KRI1', 'RSPH4A', 'PLPPR3', 'EIF5B', 'NISCH', 'HUWE1', 'TRANK1', 'PELP1', 'RBM10', 'GZF1', 'SCAF1', 'PRDM2', 'TSPOAP1', 'KAT6B', 'APLP2', 'TTBK1', 'PATL2', 'HRC', 'MYT1']))
genes.update(set(['SALL2', 'TTBK1', 'KRI1', 'ZBTB4', 'PLPPR3', 'RYR1', 'TSPOAP1', 'ZMYND15', 'RANGAP1', 'TRANK1', 'PELP1', 'RBM10', 'SCAF1', 'PRDM2', 'BBX', 'KAT6B', 'APLP2', 'HDGF', 'SMARCA2', 'MYT1']))
genes.update(set(['SP8', 'ILF3', 'EGR1', 'FUS', 'HOXB3', 'MTHFD1L', 'ONECUT2', 'PURA', 'ZFHX3', 'RFX1', 'MAP3K11', 'LURAP1L', 'CCDC6', 'HSPBP1', 'AEBP2', 'ANKRD17', 'ZNF503', 'SP3', 'TRA2A', 'DHX36']))
genes.update(set(['WT1', 'HS3ST4', 'BAG6', 'WASF2', 'KMT2E', 'SALL2', 'TRRAP', 'SF1', 'PRR12', 'FHDC1', 'RANBP9', 'WASL', 'PCDH15', 'WIPF3', 'MICALCL', 'PPP3CB', 'UBR4', 'SKOR2', 'HNRNPUL1', 'SEMA4G']))
genes.update(set(['TNRC18', 'KMT2A', 'SERTAD4', 'HDGFL2', 'MLLT6', 'PPRC1', 'SP4', 'KCNMA1', 'MLLT3', 'POU6F2', 'ERF', 'NRG3', 'ATN1', 'ZBTB4', 'SALL2', 'ZC3H3', 'BCL6B', 'LZTS2', 'NOLC1', 'GLCCI1']))
genes.update(set(['FAM76B', 'SORBS2', 'YY1', 'DLGAP3', 'CACNA1A', 'FOXG1', 'OTX1', 'SKIDA1', 'SYNGAP1', 'GATA6', 'POU4F2', 'ZIC3', 'HOXA1', 'NR4A3', 'DYRK1A', 'MAFA', 'CBX4', 'USP34', 'ONECUT1', 'MEOX2']))
genes.update(set(['POU6F2', 'ZNF384', 'MAML2', 'MAML3', 'TNS1', 'TSC22D1', 'PHLDA1', 'MEF2A', 'PHC1', 'POLG', 'GLG1', 'MAMLD1', 'HTT', 'KCNN3', 'ATXN7', 'SMARCA2', 'ATXN2', 'NCOR2']))
genes.update(set(['CNNM1', 'POU3F3', 'CHD3', 'KMT5A', 'RBM47', 'TMEM158', 'MAZ', 'MAP3K4', 'RUNX2', 'EN1', 'RBM4', 'ZIC2', 'IRF2BPL', 'ZFX', 'FOXB2', 'OLIG2', 'ZFP91', 'NEXN-AS1']))
genes.update(set(['TTBK1', 'TRANK1', 'PLPPR3', 'RYR1', 'SALL2', 'ZMYND15', 'RANGAP1', 'PELP1', 'DIEXF', 'RBM10', 'SCAF1', 'PRDM2', 'TSPOAP1', 'VIRMA', 'KAT6B', 'APLP2', 'HDGF', 'MYT1']))
genes.update(set(['DBX2', 'HMX1', 'GBX1', 'GSX1', 'HOXD8', 'HOXB9', 'CDX2', 'HOXD1', 'HOXA11', 'NKX1-1', 'HOXA4', 'DBX1', 'NKX3-2', 'HOXA1', 'MSX1', 'HOXC12', 'NKX2-6', 'DLX3']))
genes.update(set(['SHROOM4', 'DAXX', 'CHIC1', 'CERS3', 'TTBK1', 'CHGA', 'RGL2', 'ALMS1', 'FAM212A', 'CELSR2', 'IRX6', 'PRKCSH', 'CUL9', 'AEBP2', 'MAP3K9', 'HRC', 'MYT1']))
genes.update(set(['ZNF445', 'ZSCAN5DP', 'ZNF223', 'ZNF672', 'GLI4', 'ZNF213', 'ZNF174', 'ZNF75D', 'ZNF487', 'ZKSCAN5', 'ZNF781', 'ZNF350', 'ZFP69', 'ZNF853', 'ZNF552', 'PRDM9', 'ZNF629']))
genes.update(set(['ZNF215', 'ZNF536', 'ZNF554', 'ZNF319', 'ZSCAN29', 'ZNF233', 'ZNF449', 'ZNF597', 'ZNF75D', 'ZNF781', 'ZNF18', 'ZNF768', 'ZNF236', 'KLF8', 'ZNF394', 'ZNF429']))
genes.update(set(['SP8', 'TCF7L1', 'ZNF384', 'ARID1B', 'EGR1', 'RFX1', 'ZFHX3', 'EMC10', 'YEATS2', 'PURA', 'PITPNM2', 'ILF3', 'HOXB3', 'DACH1', 'WIPF1', 'HSPBP1']))
genes.update(set(['TMPRSS9', 'PRSS44', 'CTRL', 'TMPRSS6', 'TPSAB1', 'TMPRSS3', 'TMPRSS5', 'OVCH1', 'PLAT', 'PLG', 'KLK8', 'PRSS36', 'F12', 'TMPRSS13', 'TMPRSS11E', 'LPA']))
genes.update(set(['SIK2', 'FIGN', 'FBXO11', 'DAB2IP', 'MEX3C', 'SETD1A', 'ONECUT3', 'ZNF341', 'HTT', 'CCDC6', 'ZNF462', 'PCLO', 'FOXB2', 'SF3B2', 'ADGRB1', 'PEG10']))
genes.update(set(['SLC4A5', 'TRIO', 'RAD23B', 'SLITRK3', 'RALY', 'ARID1B', 'SRP68', 'TCF7L1', 'LOR', 'PURA', 'NKX1-2', 'HOXB3', 'DACH1', 'ANKS1A', 'ZSWIM8', 'KIF3C']))
genes.update(set(['LTBP3', 'BRI3BP', 'LAMP1', 'NOTCH4', 'LRP8', 'CSF1R', 'OTOF', 'LRCOL1', 'SHBG', 'ADAM33', 'HGFAC', 'BTN2A2', 'SEMA4B', 'STRCP1', 'PLXNA1', 'DPEP3']))
genes.update(set(['KLHL11', 'CHD3', 'SOX1', 'ARID1B', 'KMT5A', 'HOXA13', 'LHX2', 'CPEB3', 'EN1', 'PRPF40A', 'ZBTB8B', 'FOXB2', 'ALKBH5', 'PRKD1', 'SKOR1']))
genes.update(set(['FHOD1', 'FASLG', 'SHTN1', 'VCPIP1', 'PDCD7', 'YLPM1', 'PRR12', 'ZNF341', 'PRR11', 'PRIMA1', 'SETBP1', 'KDM6B', 'WIPF3', 'ANTXRL', 'PDE4D']))
genes.update(set(['ZBTB39', 'KLF8', 'ZNF410', 'ZNF8', 'ZNF692', 'ZSCAN29', 'ZNF597', 'ZNF324B', 'ZNF785', 'ZNF787', 'ZKSCAN4', 'ZNF275', 'REPIN1', 'ZNF137P', 'ZNF429']))
genes.update(set(['FOXP2', 'ARID3B', 'NCOA6', 'GIGYF2', 'KMT2D', 'MN1', 'POU6F2', 'GLG1', 'MAMLD1', 'AR', 'KCNN3', 'MED15', 'TOX3', 'SMARCA2', 'CREBBP']))
genes.update(set(['TMPRSS9', 'PRSS44', 'CTRL', 'TMPRSS6', 'TMPRSS3', 'TMPRSS5', 'OVCH1', 'PLAT', 'PLG', 'KLK8', 'PRSS36', 'F12', 'TMPRSS13', 'TPSAB1', 'LPA']))
genes.update(set(['SLC4A5', 'SLITRK3', 'RALY', 'ZNF746', 'ARID1B', 'PURA', 'ZFHX3', 'TCF7L1', 'LOR', 'DACH1', 'NBEA', 'CNOT3', 'ANKS1A', 'SRP68', 'KIF3C']))
genes.update(set(['BHMG1', 'MAGEA10', 'EPHB6', 'TPRXL', 'SALL1', 'CHD8', 'ZNF865', 'REXO1', 'CAMTA2', 'SHE', 'SYN2', 'KLHL1', 'NR1D1', 'RANBP10', 'NKX6-1']))
genes.update(set(['POU3F3', 'WDR26', 'CAPNS1', 'HOXD9', 'FEZF2', 'ARID1B', 'FUBP1', 'CACNG8', 'POU4F2', 'SHOX2', 'LRCH2', 'MAF', 'TAL1', 'SOX1', 'YEATS2']))
genes.update(set(['ZNF536', 'ZNF554', 'ZNF319', 'ZSCAN29', 'ZNF233', 'ZNF449', 'ZNF597', 'ZNF75D', 'ZNF781', 'ZNF768', 'ZNF236', 'KLF8', 'ZNF394', 'ZNF429']))
genes.update(set(['LHX4', 'ARX', 'ALX1', 'DRGX', 'UNCX', 'ESX1', 'SHOX2', 'VSX1', 'PHOX2A', 'OTP', 'PRRX1', 'MIXL1', 'RAX', 'PROP1']))
#genes.update(set(['NAF1', 'FHOD1', 'LMTK3', 'SHTN1', 'YLPM1', 'PRR12', 'ZNF341', 'PRR11', 'FASLG', 'SETBP1', 'VCPIP1', 'WIPF3', 'ANTXRL', 'PDE4D']))
#genes.update(set(['EHMT2', 'WIZ', 'VSIG10', 'SENP3', 'CELSR2', 'CDK11B', 'SLC24A1', 'SMYD5', 'ZFP91', 'FBLN2', 'DCAF8', 'SMARCA4', 'CLSTN1', 'PHF23']))
#genes.update(set(['FASLG', 'PRIMA1', 'MYPOP', 'VCPIP1', 'SETD1B', 'PRR12', 'ZNF341', 'PRR11', 'DBN1', 'SETBP1', 'KDM6B', 'WIPF3', 'BCL9L', 'PDE4D']))
#genes.update(set(['ZBTB7B', 'TRERF1', 'MYPOP', 'RAPH1', 'RNF44', 'FMN1', 'PELP1', 'SHROOM3', 'ACR', 'SETD1B', 'PRIMA1', 'KDM6B', 'FASLG', 'PRRG2']))
#genes.update(set(['DCHS1', 'PTPRN2', 'PCSK9', 'TMUB1', 'LRP5', 'CELSR2', 'ADAMTS8', 'PCDHAC2', 'CNPY3', 'ADAM33', 'LAMP1', 'SEMA4B', 'STRCP1', 'FUCA2']))
#genes.update(set(['NEUROD1', 'ARIH2', 'HINFP', 'HOMEZ', 'KIF21A', 'NAP1L2', 'ACIN1', 'ZBTB7C', 'CSRNP3', 'KDM2A', 'SEC63', 'FKBP8', 'ZBTB9', 'NPM2']))
#genes.update(set(['POU3F3', 'SOX21', 'FBXO41', 'ARX', 'ZNF358', 'ARID1A', 'MSH3', 'HOXA13', 'JUND', 'TBX2', 'POU4F1', 'PHOX2B', 'MNX1', 'LHFPL3']))
#genes.update(set(['TRIO', 'ANO8', 'SLITRK3', 'RALY', 'ARID1B', 'PURA', 'ZFHX3', 'TCF7L1', 'LOR', 'DACH1', 'NBEA', 'CNOT3', 'ANKS1A', 'SRP68']))
#genes.update(set(['KLHL11', 'CHD3', 'SOX1', 'ARID1B', 'DEAF1', 'CPEB3', 'EN1', 'FOXO1', 'MED14', 'ZBTB8B', 'FOXB2', 'ALKBH5', 'SKOR1']))
#genes.update(set(['SP8', 'FUS', 'ARID1B', 'ONECUT2', 'PURA', 'ZFHX3', 'RFX1', 'MAP3K11', 'LURAP1L', 'AEBP2', 'HNRNPAB', 'HOXB3', 'DHX36']))
#genes.update(set(['TMPRSS9', 'PRSS56', 'CTRL', 'TPSD1', 'PRSS12', 'TMPRSS3', 'PLG', 'TMPRSS13', 'ST14', 'PLAU', 'TMPRSS15', 'KLK9', 'LPA']))
#genes.update(set(['OTUD7A', 'RBM4', 'ARX', 'CDK13', 'NIPA1', 'ASCL1', 'RPL14', 'RBM24', 'RBM23', 'DACH1', 'SOX3', 'NRG3', 'TMEM158']))
#genes.update(set(['TCF7L1', 'SNX27', 'RALY', 'ARID1B', 'RFX1', 'EMC10', 'LOR', 'PURA', 'KCNN2', 'HOXB3', 'DACH1', 'YEATS2', 'KIF3C']))
#genes.update(set(['SCAF4', 'SOCS7', 'ZNF384', 'MAML3', 'FBXO11', 'MEF2A', 'PAXIP1', 'POU6F2', 'AAK1', 'HTT', 'KCNN3', 'ATXN7', 'ATXN2']))
#genes.update(set(['TMPRSS9', 'PRSS56', 'PLAU', 'F9', 'PRSS12', 'TMPRSS3', 'PLG', 'TMPRSS13', 'ST14', 'KLK5', 'TMPRSS15', 'TPSAB1', 'LPA']))
#genes.update(set(['ZNF672', 'MZF1', 'ZNF200', 'ZSCAN29', 'ZNF487', 'KLF7', 'ZNF669', 'ZKSCAN4', 'RBAK', 'ZNF775', 'ZNF137P', 'PRDM9', 'ZNF629']))
#genes.update(set(['POU3F3', 'SRSF1', 'FEZF2', 'ONECUT3', 'NPAS3', 'EVX2', 'CACNG8', 'HNRNPL', 'FOXF1', 'MAF', 'TAL1', 'ZNF503', 'SOX1']))
#genes.update(set(['ZNF777', 'ZNF628', 'ZNF350', 'ZNF316', 'KLF14', 'ZNF274', 'ZNF853', 'ZNF768', 'ZNF500', 'ZNF416', 'ZSCAN5A', 'ZFP69B', 'ZNF629']))
#genes.update(set(['SIK3', 'MAML3', 'FAM157B', 'TBP', 'ZFHX3', 'FAM155A', 'POLG', 'MED12', 'EP400', 'TMEM229A', 'DENND4B', 'ATXN7', 'THAP11']))
#genes.update(set(['DBX2', 'DBX1', 'HOXC12', 'HOXC4', 'HOXA11', 'HOXD1', 'HOXA6', 'HOXB9', 'HOXA3', 'HOXA1', 'PDX1', 'MNX1']))
#genes.update(set(['SALL2', 'WASF2', 'FMNL1', 'ZC3H4', 'HOXB4', 'SF3B2', 'ESRP2', 'HTT', 'WASL', 'PCDH15', 'PPP3CB', 'SF1']))
#genes.update(set(['SRRM2', 'SERTAD4', 'TPRXL', 'ZC3H3', 'RNF126', 'ERF', 'CUL4B', 'ZBTB4', 'KLHL1', 'RANBP10', 'RLIM', 'NKX6-1']))
#genes.update(set(['HMX3', 'TMEM64', 'TBL1X', 'RBM47', 'FOXD3', 'SP5', 'EVX2', 'CKAP4', 'ZFP36L2', 'PITX3', 'FOXC1', 'NKX6-1']))
#genes.update(set(['SP8', 'CDK13', 'HOXA13', 'CPEB3', 'EVX2', 'TMEM200C', 'GATA6', 'PITX3', 'SOX1', 'IRX1', 'CKAP4', 'NKX6-1']))
#genes.update(set(['KLF12', 'ZNF217', 'ZNF202', 'MYNN', 'ZNF826P', 'ZNF137P', 'ZNF48', 'KLF5', 'ZNF433', 'ZNF287', 'ZNF76', 'ZNF396']))
#genes.update(set(['POU3F3', 'WDR26', 'CAPNS1', 'HOXD9', 'RALY', 'ARID1B', 'FOXD1', 'CACNG8', 'LRCH2', 'TAL1', 'SOX1', 'YEATS2']))
#genes.update(set(['ZNF837', 'ZNF775', 'ZNF542P', 'ZNF629', 'MZF1', 'ZNF648', 'ZSCAN9', 'ZFP69', 'ZNF671', 'ZNF37A', 'ZNF282', 'ZNF467']))
#genes.update(set(['POU3F3', 'IFFO2', 'ADCY1', 'CCM2L', 'ARID1B', 'NPAS3', 'EVX2', 'CACNG8', 'ONECUT3', 'HNRNPL', 'FOXD1', 'CCDC85C']))
#genes.update(set(['ELN', 'BHLHE41', 'FOXD3', 'MAP3K4', 'CCDC177', 'EN1', 'ZIC2', 'FOXF2', 'SKIDA1', 'PHLPP1', 'OLIG2', 'GDF11']))
#genes.update(set(['SP8', 'FUS', 'RALY', 'ARID1B', 'AEBP2', 'ONECUT2', 'RFX1', 'ZFHX3', 'SYNGAP1', 'LURAP1L', 'HOXB3', 'ZNF281']))
#genes.update(set(['VCPIP1', 'TRERF1', 'MYPOP', 'KDM6B', 'FMN1', 'RAPH1', 'SETD1B', 'PRIMA1', 'SETBP1', 'PWWP2A', 'FASLG', 'PDE4D']))
#genes.update(set(['NEUROD2', 'YTHDC1', 'PRM3', 'SKIDA1', 'CENPB', 'TEX2', 'SPRYD3', 'PRDM2', 'DCAF8L2', 'KDM2A', 'FKBP8', 'HRC']))
#genes.update(set(['MIDN', 'GBX1', 'PRDM12', 'ATP6AP1', 'TGFBR1', 'TMEM200C', 'SOWAHA', 'ZIC5', 'GATA4', 'IRX1', 'NBEAL2', 'NEXN-AS1']))
#genes.update(set(['ZNF445', 'ZNF334', 'ZNF629', 'ZSCAN5DP', 'ZNF749', 'ZNF24', 'ZNF781', 'ZKSCAN2', 'ZNF853', 'GLI4', 'ZNF473', 'PRDM9']))
#genes.update(set(['KLF12', 'ZNF217', 'ZNF195', 'ZNF202', 'MYNN', 'ZNF826P', 'ZNF48', 'KLF5', 'ZNF433', 'ZNF287', 'ZNF76', 'ZNF396']))
#genes.update(set(['POU3F3', 'TWIST1', 'MLLT6', 'GSX1', 'NPAS3', 'HOXA10', 'ZSWIM6', 'ZIC2', 'LRCH2', 'MAFA', 'FOXD1']))
#genes.update(set(['DDX17', 'FHDC1', 'KMT2E', 'SCAF11', 'ZFHX4', 'ZNF462', 'DMRTB1', 'RNF215', 'WAS', 'PEG10', 'LMOD2']))
#genes.update(set(['FOXJ2', 'ZNF384', 'MAML3', 'SOCS7', 'IER5L', 'POU6F2', 'MEF2A', 'HTT', 'KCNN3', 'ATXN7', 'ATXN2']))
#genes.update(set(['RALY', 'ARID1B', 'RFX1', 'ZFHX3', 'TCF7L1', 'PURA', 'DACH1', 'ANKS1A', 'SRP68', 'ZNF281', 'COL17A1']))
#genes.update(set(['KMT2A', 'EPHB6', 'ZBTB4', 'SYN2', 'NFE2L1', 'POU6F2', 'ATN1', 'SHE', 'NRG3', 'LZTS2', 'CLASRP']))
#genes.update(set(['USF3', 'CACNA1A', 'ASCL1', 'ZFHX3', 'RAI1', 'MAMLD1', 'MED15', 'SATB1', 'TBP', 'ST6GALNAC5', 'NUMBL']))
#genes.update(set(['CHD1', 'FOXP2', 'MAML3', 'KMT2D', 'ARID1B', 'BMP2K', 'MED12', 'ATN1', 'VEZF1', 'AMOT', 'ATXN1']))
#genes.update(set(['ZBTB16', 'ZNF619', 'ZNF302', 'ZNF8', 'ZNF24', 'ZNF552', 'GLI4', 'ZNF577', 'ZNF114', 'ZNF439', 'ZIM3']))
#genes.update(set(['BAG6', 'KMT2E', 'FNBP4', 'TPRN', 'ENAH', 'DAB2IP', 'USP51', 'PRR12', 'WASH2P', 'SCAF1', 'APBB1IP']))
#genes.update(set(['LHX4', 'ARX', 'ALX1', 'DRGX', 'UNCX', 'ESX1', 'SHOX2', 'PHOX2A', 'OTP', 'PRRX1', 'RAX']))
#genes.update(set(['OTUD7A', 'KDM1A', 'GBX1', 'ARX', 'ZFR', 'DMRTA1', 'COPS6', 'EOMES', 'FBRS', 'DACH1', 'NKX2-4']))
#genes.update(set(['CHD1', 'HELLS', 'CHD3', 'SRCAP', 'RAD54B', 'CHD6', 'SHPRH', 'INO80', 'SMARCAD1', 'SMARCA1', 'SMARCA4']))
#genes.update(set(['PRSS44', 'TMPRSS6', 'KLK12', 'TMPRSS5', 'OVCH1', 'PLAT', 'HPN', 'F12', 'PRSS8', 'F10', 'TMPRSS11E']))
#genes.update(set(['KLHL11', 'CHD3', 'DEAF1', 'ZFHX3', 'EN1', 'FOXO1', 'MED14', 'ZBTB8B', 'FOXB2', 'ALKBH5', 'SKOR1']))
#genes.update(set(['ZNF597', 'KLF8', 'ZNF785', 'ZSCAN29', 'ZNF8', 'ZBTB39', 'ZNF787', 'ZKSCAN4', 'ZNF275', 'REPIN1', 'ZNF429']))
#genes.update(set(['ZNF628', 'ZNF350', 'ZNF316', 'ZFP69B', 'ZNF853', 'ZNF768', 'ZNF274', 'ZNF500', 'ZSCAN5A', 'ZNF671', 'ZNF629']))
#genes.update(set(['CHD1', 'HELLS', 'CHD3', 'SRCAP', 'RAD54B', 'CHD6', 'SHPRH', 'INO80', 'SMARCAD1', 'SMARCA1', 'SMARCA4']))
#genes.update(set(['AKNA', 'FAM9B', 'SKIDA1', 'KCNH8', 'PPARGC1B', 'CENPB', 'PRM3', 'KAT6B', 'RSPH6A', 'RANGAP1', 'NPM2']))
#genes.update(set(['ZNF445', 'ZNF324B', 'ZNF619', 'ZNF233', 'ZNF697', 'ZFP41', 'ZKSCAN5', 'ZNF513', 'RBAK', 'ZNF500', 'ZSCAN10']))
#genes.update(set(['ANP32E', 'PHACTR1', 'ZNF316', 'RTN4', 'HMGB2', 'SENP3', 'FBXO3', 'ABCC9', 'PIAS4', 'KANK1', 'ARID4B']))
#genes.update(set(['ZBTB16', 'ZNF619', 'ZNF302', 'ZNF8', 'ZNF487', 'ZNF24', 'ZNF552', 'GLI4', 'ZNF577', 'ZNF439']))
#genes.update(set(['SHROOM4', 'PLEKHG5', 'TPRN', 'RGL2', 'SLC24A1', 'CNPY4', 'P3H3', 'ZBTB47', 'TSPOAP1', 'MYO15B']))
#genes.update(set(['ZNF445', 'ZNF202', 'MZF1', 'ZSCAN18', 'ZBED9', 'ZNF449', 'ZNF496', 'ZKSCAN5', 'ZNF274', 'ZSCAN10']))
#genes.update(set(['BHMG1', 'CHD9', 'JCAD', 'SETD1A', 'DCAF5', 'SNAPC4', 'BRDT', 'BCL6B', 'NRG2', 'TMEM40']))
#genes.update(set(['KLF11', 'ZNF202', 'GLI4', 'ZNF217', 'KLF5', 'HKR1', 'KLF1', 'ZBTB48', 'ZBTB24', 'ZNF560']))
#genes.update(set(['CHD5', 'MLLT6', 'MAZ', 'EN1', 'USP11', 'ZMIZ2', 'KCNA4', 'ALKBH5', 'ASXL3', 'STUM']))
#genes.update(set(['ZNF202', 'ZNF200', 'ZNF316', 'ZNF146', 'ZNF75A', 'ZNF487', 'KLF7', 'ZNF18', 'ZSCAN20', 'PRDM9']))
#genes.update(set(['ZFHX3', 'MAP3K9', 'CCND1', 'TIMELESS', 'ADRA2B', 'SCAF1', 'WDR70', 'ARID4B', 'C11orf95', 'DHX37']))
#genes.update(set(['ARX', 'ALX1', 'DRGX', 'UNCX', 'ESX1', 'SHOX2', 'OTP', 'PHOX2A', 'PRRX1', 'RAX']))
#genes.update(set(['DBX2', 'HOXB9', 'HOXD8', 'HOXD1', 'HOXA11', 'HOXC12', 'HOXA1', 'HOXA4', 'MNX1', 'DBX1']))
#genes.update(set(['SP8', 'TGFBR1', 'GSX2', 'BHLHE22', 'FOXD3', 'SP5', 'EVX2', 'ZIC2', 'PHLPP1', 'PRKD1']))
#genes.update(set(['ICK', 'RPS6KA5', 'BRSK1', 'CAMK1', 'CAMK2B', 'MELK', 'PRKACG', 'DCLK3', 'MAPKAPK5', 'PSKH1']))
#genes.update(set(['B4GALNT3', 'SKIDA1', 'TRIM26', 'AKNA', 'RANGAP1', 'CENPB', 'SETD1A', 'MICAL3', 'USP34', 'NPM2']))
#genes.update(set(['KLF12', 'ZNF613', 'ZNF205', 'ZNF672', 'ZNF316', 'ZNF826P', 'ZBTB47', 'ZNF500', 'ZNF799', 'ZNF629']))
#genes.update(set(['PDGFRB', 'CSK', 'EPHA8', 'FRK', 'EGFR', 'ZAP70', 'JAK1', 'TIE1', 'FGFR2', 'FLT3']))
#genes.update(set(['KLHL11', 'CHD3', 'DEAF1', 'ZFHX3', 'EN1', 'FOXO1', 'ZBTB8B', 'FOXB2', 'ALKBH5', 'SKOR1']))
#genes.update(set(['SP9', 'DMRTA1', 'IRS4', 'SLC12A2', 'ASXL2', 'HOXD13', 'CASKIN1', 'C6orf223', 'PHOX2B', 'ZFP91']))
#genes.update(set(['SOCS7', 'MAML3', 'NONO', 'FBXO11', 'SF1', 'MEF2A', 'HTT', 'KCNN3', 'ATXN7', 'PHLPP1']))
#genes.update(set(['SOX21', 'DMRTA2', 'ARID1B', 'DEAF1', 'JUND', 'REPS2', 'C6orf223', 'NRG3', 'PHOX2B', 'MNX1']))
#genes.update(set(['EDA', 'COL18A1', 'COLQ', 'ERFE', 'COL4A1', 'COL16A1', 'EMID1', 'COL20A1', 'COL12A1', 'COL17A1']))
#genes.update(set(['SKIDA1', 'FAM9A', 'KCNH8', 'PPARGC1B', 'CENPB', 'ANP32C', 'PRM3', 'KAT6B', 'RANGAP1', 'NPM2']))
#genes.update(set(['PTPN21', 'PTPRQ', 'PTPRG', 'PTPRD', 'PTPRC', 'PTPRB', 'PTPRA', 'PTPRO', 'PTPRJ', 'PTPRH']))
#genes.update(set(['KIF16B', 'STARD9', 'KIF6', 'KIF7', 'KIF14', 'KIF15', 'KIF13B', 'KIF17', 'KIF2B', 'KIFC1']))
#genes.update(set(['RLIM', 'NFE2L1', 'BHLHE22', 'NBEA', 'EDC4', 'SNAPC4', 'DACH1', 'TSHZ3', 'MEX3B', 'TNRC18']))

all_prots = list(genes)



for p in all_prots:
    print p

# dilute similar proteins
diluted_list = []
pb = Progressbar('Diluting proteins...')
i = 0
for protein_name1 in all_prots:
    i += 1
    pb.update_progress(i, len(all_prots))
    redundantProt = False
    for protein_name in diluted_list:
        #print '%s: Checking similarity of %s and %s' % (kmer, protein_name, prot.geneName)
        if not utils.proteins_are_dissimilar(protein_name, protein_name1,
                                             protein_seq[protein_name], protein_seq[protein_name1]):
            redundantProt = True
            break
    if not redundantProt:
        diluted_list.append(protein_name1)
    else:
        print 'Ignoring redundant protein: %s' % protein_name1
    redundantProt = False

for p in diluted_list:
    print p

