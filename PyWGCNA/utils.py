import pickle
import os

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
