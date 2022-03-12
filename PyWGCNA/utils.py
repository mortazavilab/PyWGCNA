import pickle
import os
import biomart

from PyWGCNA.comparison import *

# bcolors
HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKCYAN = '\033[96m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'


# read WGCNA obj
def readWGCNA(file):
    """
    Read a WGCNA from a saved pickle file.

    :param file: Name / path of WGCNA object
    :type file: str

    :return: PyWGCNA object
    :rtype: PyWGCNA class
    """
    if not os.path.isfile(file):
        raise ValueError('WGCNA object not found at given path!')

    picklefile = open(file, 'rb')
    wgcna = pickle.load(picklefile)

    print(f"{BOLD}{OKBLUE}Reading WGCNA done!{ENDC}")
    return wgcna


# compare two WGCNA
def compareWGCNA(WGCNA1, WGCNA2):
    """
    Compare two WGCNAs
                
    :param WGCNA1: first WGCNA object
    :type WGCNA1: PyWGCNA class
    :param WGCNA2: second WGCNA object
    :type WGCNA2: PyWGCNA class

    :return: compare object
    :rtype: Compare class
    """
    compare = Comparison(name1=WGCNA1.name, name2=WGCNA2.name,
                         geneModule1=WGCNA1.geneModules, geneModule2=WGCNA2.geneModules)
    compare.compareWGCNA()

    return compare


# compare WGCNA to single cell
def compareSingleCell(WGCNA, sc):
    """
    Compare WGCNA and gene marker from single cell experiment

    :param WGCNA: WGCNA object
    :type WGCNA: PyWGCNA class
    :param sc: gene marker table which has ....
    :type sc: pandas dataframe

    :return: compare object
    :rtype: Compare class

    """
    compare = Comparison(name1=WGCNA.name, geneModule1=WGCNA.geneModules,
                         geneMarker=sc, sc=True)
    compare.compareSingleCell()

    return compare


def getGeneList(dataset='mmusculus_gene_ensembl', attributes=None):
    """
    get dictionary that map gene ensembl id to gene name from biomart
    
    
    :param dataset: name of the dataset we used from biomart; mouse: mmusculus_gene_ensembl and human: hsapiens_gene_ensembl
        you can find more information here: https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#selecting-a-biomart-database-and-dataset
    :type dataset: string
    :param attributes: List the types of data we want
    :type attributes: list
    
    :return: table extracted from biomart related to the datasets including information from attributes
    :rtype: pandas dataframe
    """
    # Set up connection to server
    if dataset == 'mmusculus_gene_ensembl':
        attributes = ['ensembl_transcript_id', 'mgi_symbol', 'ensembl_gene_id']
    if dataset == 'hsapiens_gene_ensembl':
        attributes = ['ensembl_transcript_id', 'hgnc_symbol', 'ensembl_gene_id']

    server = biomart.BiomartServer('http://uswest.ensembl.org/biomart')
    mart = server.datasets[dataset]

    # Get the mapping between the attributes
    response = mart.search({'attributes': attributes})
    data = response.raw.data.decode('ascii')

    geneInfo = pd.DataFrame(columns=attributes)
    # Store the data in a dict
    for line in data.splitlines():
        line = line.split('\t')
        dict = {attributes[0]: line[0],
                attributes[1]: line[1],
                attributes[2]: line[2]}
        geneInfo = geneInfo.append(dict, ignore_index=True)

    return geneInfo


# read comparison obj
def readComparison(file):
    """
    Read a comparison from a saved pickle file.

    :param file: Name / path of comparison object
    :type file: string

    :return: comparison object
    :rtype: comparison class
    """
    if not os.path.isfile(file):
        raise ValueError('Comparison object not found at given path!')

    picklefile = open(file, 'rb')
    comparison = pickle.load(picklefile)

    print(f"{BOLD}{OKBLUE}Reading comparison done!{ENDC}")
    return comparison
