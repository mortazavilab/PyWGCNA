import numpy as np
import pandas as pd
import os
import anndata as ad


# remove runtime warning (divided by zero)
np.seterr(divide='ignore', invalid='ignore')


class GeneExp:
    """
    A class used to creat gene expression anndata along data trait including both genes and samples information.

    :param species: species of the data you use i.e mouse, human
    :type species: str
    :param level: which type of data you use including gene, transcript (default: gene)
    :type level: str
    :param anndata: if the expression data is in anndata format you should pass it through this parameter. X should be expression matrix. var is a sample information and obs is a gene information.
    :param anndata: anndata
    :param geneExpr: expression matrix which genes are in the rows and samples are columns
    :type geneExpr: pandas dataframe
    :param geneExpPath: path of expression matrix
    :type geneExpPath: str
    :param sep: separation symbol to use for reading data in geneExpPath properly
    :type sep: str
    """

    def __init__(self, species=None, level='gene',
                 anndata=None, geneExp=None,
                 geneExpPath=None, sep=' '):
        self.species = species
        self.level = level
        if geneExpPath is not None:
            if not os.path.isfile(geneExpPath):
                raise ValueError("file does not exist!")
            else:
                expressionList = pd.read_csv(geneExpPath, sep=sep)
        elif geneExp is not None:
            if isinstance(geneExp, pd.DataFrame):
                expressionList = geneExp
            else:
                raise ValueError("geneExp is not data frame!")
        elif anndata is not None:
            if isinstance(anndata, ad.AnnData):
                self.geneExpr = anndata
                return
            else:
                raise ValueError("geneExp is not data frame!")
        else:
            raise ValueError("all type of input can not be empty at the same time!")

        column = 'id'
        if level == 'gene':
            column = 'gene_id'
        elif level == 'transcript':
            column = 'transcript_id'

        geneInfo = pd.DataFrame(expressionList.values[:, 0], columns=[column],
                                index=expressionList.iloc[:, 0])

        sampleInfo = pd.DataFrame(range(len(expressionList.columns[1:])), columns=['sample_id'],
                                  index=expressionList.columns[1:])

        expressionList.index = expressionList.iloc[:, 0]  # gene_id
        # drop gene id columns
        expressionList = expressionList.drop([expressionList.columns[0]], axis=1)

        self.geneExpr = ad.AnnData(X=expressionList, obs=geneInfo, var=sampleInfo)

    @staticmethod
    def updateGeneInfo(expr, geneInfo=None, path=None, sep=' ', order=True, level='gene'):
        """
        add/update genes info in expr anndata

        :param expr: expression data
        :type expr: anndata
        :param geneInfo: gene information table you want to add to your data
        :type geneInfo: pandas dataframe
        :param path: path of geneInfo
        :type path: str
        :param sep: separation symbol to use for reading data in path properly
        :type sep: str
        :param order: if you want to update/add gene information by keeping the order as the same as data. if you want to add gene infor from biomart you should set this to be false. (default: TRUE)
        :type order: bool
        :param level: indicated the expression data is at gene level or transcript level
        :type level: str
        """
        if path is not None:
            if not os.path.isfile(path):
                raise ValueError("path does not exist!")
            geneInfo = pd.read_csv(path, sep=sep)
        elif geneInfo is not None:
            if not isinstance(geneInfo, pd.DataFrame):
                raise ValueError("geneInfo is not pandas dataframe!")
        else:
            raise ValueError("path and geneInfo can not be empty at the same time!")

        if order:
            geneInfo.index = expr.obs
            expr.obs = pd.concat([geneInfo, expr.obs], axis=1)
            expr.obs = expr.obs.loc[:, ~expr.obs.columns.duplicated()]
        else:
            name = 'ensembl_gene_id'
            replace = 'gene_id'
            if level == 'transcript':
                name = 'ensembl_transcript_id'
                replace = 'transcript_id'
            if 'external_gene_name' in geneInfo.columns:
                geneInfo.rename(columns={'external_gene_name': 'gene_name', name: replace}, inplace=True)
            else:
                geneInfo.rename(columns={name: replace}, inplace=True)
            expr.obs.gene_id = expr.obs.gene_id.str.split('\\.', expand=True)[0]
            expr.obs.index.name = None
            rmv = [x for x in geneInfo.columns if x in expr.obs.columns]
            rmv.remove(replace)
            expr.obs.drop(rmv, axis=1, inplace=True)
            expr.obs = expr.obs.merge(geneInfo, on=replace, how='left')
            expr.obs.index = expr.obs[replace]
        return expr

    @staticmethod
    def updateMetadata(expr, metaData=None, path=None, sep=' '):
        """
        add/update metadata in expr anndata

        :param expr: expression data
        :type expr: anndata
        :param metaData: Sample information table you want to add to your data
        :type metaData: pandas dataframe
        :param path: path of metaData
        :type path: str
        :param sep: separation symbol to use for reading data in path properly
        :type sep: str
        """
        if path is not None:
            if not os.path.isfile(path):
                raise ValueError("path does not exist!")
            metaData = pd.read_csv(path, sep=sep)
        elif metaData is not None:
            if not isinstance(metaData, pd.DataFrame):
                raise ValueError("meta data is not pandas dataframe!")
        else:
            raise ValueError("path and metaData can not be empty at the same time!")

        metaData.index = expr.var.index
        expr.var = pd.concat([metaData, expr.var], axis=1)
        expr.var = expr.var.loc[:, ~expr.var.columns.duplicated()]

        return expr
