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
    :param sep: str
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

        geneInfo = pd.DataFrame(expressionList.values[:, 0], columns=['gene_id'],
                                index=expressionList.iloc[:, 0])

        sampleInfo = pd.DataFrame(range(len(expressionList.columns[1:])), columns=['sample_id'],
                                  index=expressionList.columns[1:])

        expressionList.index = expressionList.iloc[:, 0]  # gene_id
        # drop gene id columns
        expressionList = expressionList.drop([expressionList.columns[0]], axis=1)

        self.geneExpr = ad.AnnData(X=expressionList, obs=geneInfo, var=sampleInfo)

    def updateGeneInfo(self, geneInfo=None, path=None, sep=' ', order=True):
        """
        add/update genes info in geneExp anndata

        :param geneInfo: gene information table you want to add to your data
        :type geneInfo: pandas dataframe
        :param path: path of geneInfo
        :type path: str
        :param sep: separation symbol to use for reading data in path properly
        :type sep: str
        :param order: if you want to update/add gene information by keeping the order as the same as data. if you want to add gene infor from biomart you should set this to be false. (default: TRUE)
        :type order: bool
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
            geneInfo.index = self.geneExpr.obs
            self.geneExpr.obs = pd.concat([geneInfo, self.geneExpr.obs], axis=1)
            self.geneExpr.obs = self.geneExpr.obs.loc[:, ~self.geneExpr.obs.columns.duplicated()]
        else:
            name = 'ensembl_gene_id'
            if self.level == 'transcript':
                name = 'ensembl_transcript_id'
            columns = geneInfo.columns
            columns.remove(name)
            for column in columns:
                if 'symbol' in column:
                    column = 'gene_name'
                self.geneExpr.obs[column] = ""
            for i in range(self.geneExpr.obs.shape[0]):
                if self.geneExpr.obs['gene_id'][i] in self.geneExpr.obs[name]:
                    for column in columns:
                        if 'symbol' in column:
                            column = 'gene_name'
                        self.geneExpr.obs[column][i] = geneInfo[
                            geneInfo[name] == self.geneExpr.obs['gene_id'][i], column]

    def updateMetadata(self, metaData, path, sep=' '):
        """
        add/update metadata in geneExp anndata

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

        metaData.index = self.geneExpr.var.index
        self.geneExpr.var = pd.concat([metaData, self.geneExpr.var], axis=1)
        self.geneExpr.var = self.geneExpr.var.loc[:, ~self.geneExpr.var.columns.duplicated()]
