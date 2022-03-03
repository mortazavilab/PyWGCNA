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
                
    Parameters
    ----------
    file (str): Name / path of WGCNA object

    Returns
    -------
    WGCNA (PyWGCNA): WGCNA object

    """""
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
                
    Parameters
    ----------
    WGCNA1 (WGCNA class): first WGCNA object
    WGCNA2 (WGCNA class): second WGCNA object

    Returns
    -------
    compare (Compare class)

    """""
    compare = Comparison(name1=WGCNA1.name, name2=WGCNA2.name,
                         geneModule1=WGCNA1.geneModules, geneModule2=WGCNA2.geneModules)
    compare.compareWGCNA()

    return compare


# compare WGCNA to single cell
def compareSingleCell(WGCNA, sc):
    """
        Compare WGCNA and gene marker from single cell experiment

        Parameters
        ----------
        WGCNA (WGCNA class): WGCNA object
        sc (dataframe): gene marker table which has ....

        Returns
        -------
        compare (Compare class)

        """""
    compare = Comparison(name1=WGCNA.name, geneModule1=WGCNA.geneModules,
                         geneMarker=sc, sc=True)
    compare.compareSingleCell()

    return compare


def getGeneList(dataset='mmusculus_gene_ensembl', attributes=None):
    """
    get dictionary that map gene ensembl id to gene name from biomart
    Parameters
    ----------
    dataset: name of the dataset we used from biomart; mouse: mmusculus_gene_ensembl and human: hsapiens_gene_ensembl
    you can find more information here: https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#selecting-a-biomart-database-and-dataset
    
    attributes: List the types of data we want
    
    Returns
    -------
    dictionary contain gene id as a key and gene name as a value
    """""
    # Set up connection to server
    if dataset == 'mmusculus_gene_ensembl':
        attributes = ['ensembl_transcript_id', 'mgi_symbol', 'ensembl_gene_id', 'ensembl_peptide_id']
    if dataset == 'hsapiens_gene_ensembl':
        attributes = ['ensembl_transcript_id', 'hgnc_symbol', 'ensembl_gene_id', 'ensembl_peptide_id']

    server = biomart.BiomartServer('http://uswest.ensembl.org/biomart')
    mart = server.datasets[dataset]

    # Get the mapping between the attributes
    response = mart.search({'attributes': attributes})
    data = response.raw.data.decode('ascii')

    ensembl_to_genesymbol = {}
    # Store the data in a dict
    for line in data.splitlines():
        line = line.split('\t')
        # The entries are in the same order as in the `attributes` variable
        transcript_id = line[0]
        gene_symbol = line[1]
        ensembl_gene = line[2]
        ensembl_peptide = line[3]

        # Some of these keys may be an empty string. If you want, you can
        # avoid having a '' key in your dict by ensuring the
        # transcript/gene/peptide ids have a nonzero length before
        # adding them to the dict
        ensembl_to_genesymbol[transcript_id] = gene_symbol
        ensembl_to_genesymbol[ensembl_gene] = gene_symbol
        ensembl_to_genesymbol[ensembl_peptide] = gene_symbol

    return ensembl_to_genesymbol


# read comparison obj
def readComparison(file):
    """
    Read a comparison from a saved pickle file.

    Parameters
    ----------
    file (str): Name / path of comparison object

    Returns
    -------
    comparison: comparison object

    """""
    if not os.path.isfile(file):
        raise ValueError('Comparison object not found at given path!')

    picklefile = open(file, 'rb')
    comparison = pickle.load(picklefile)

    print(f"{BOLD}{OKBLUE}Reading comparison done!{ENDC}")
    return comparison
