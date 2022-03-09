import numpy as np
import pandas as pd
import os
import anndata as ad


# remove runtime warning (divided by zero)
np.seterr(divide='ignore', invalid='ignore')


class GeneExp:
    """
        A class used to creat gene expression matrix along data trait including both genes and samples information.

        Attributes
        ----------
        geneExpr: anndata
            gene expression matrix along sample information in the var (i.e. age, sex, tissue, genotype)
             and gene information in the obs(i.e. gene name, gene id, gene type)

        Methods
        -------
        updateMetadata(path)
            add/update genes info
        addSampleInfo(path)
            add/update samples info
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

        geneInfo = pd.DataFrame(expressionList.iloc[:, 0], columns=['gene_id'])

        sampleInfo = pd.DataFrame(range(len(expressionList.columns[1:])), columns=['sample_id'], index=expressionList.columns[1:])

        expressionList.index = expressionList.iloc[:, 0]  # gene_id
        # drop gene id columns
        expressionList = expressionList.drop([expressionList.columns[0]], axis=1)

        self.geneExpr = ad.AnnData(X=expressionList, obs=geneInfo, var=sampleInfo)

    def updateGeneInfo(self, geneInfo=None, path=None, sep=' ', order=True):
        if path is not None:
            if not os.path.isfile(path):
                raise ValueError("file does not exist!")
            geneInfo = pd.read_csv(path, sep=sep)

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

    def updateMetadata(self, path, sep=' '):
        if not os.path.isfile(path):
            raise ValueError("file does not exist!")

        samples = pd.read_csv(path, sep=sep)
        samples.index = self.geneExpr.var.index
        self.geneExpr.var = pd.concat([samples, self.geneExpr.var], axis=1)
        self.geneExpr.var = self.geneExpr.var.loc[:, ~self.geneExpr.var.columns.duplicated()]
