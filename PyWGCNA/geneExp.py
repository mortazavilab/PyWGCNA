import numpy as np
import pandas as pd
import os


# remove runtime warning (divided by zero)
np.seterr(divide='ignore', invalid='ignore')


class GeneExp:
    """
        A class used to creat gene expression matrix along data trait including both genes and samples information.

        Attributes
        ----------
        geneInfo : Dataframe
            keep information about genes (i.e. gene name, gene id, gene type)

        sampleInfo : Dataframe
            keep information about samples (i.e. age, sex, tissue, genotype)

        geneExpression : Dataframe
            gene expression matrix which columns are genes and rows are samples

        Methods
        -------
        addGeneInfo(path)
            add genes info
        addSampleInfo(path)
            add samples info
        getGeneExp()
            returns gene expression matrix
        """

    def __init__(self, geneExpPath, sep=' '):
        if not os.path.isfile(geneExpPath):
            raise ValueError("file does not exist!")

        self.expressionList = pd.read_csv(geneExpPath, sep=sep)

        self.geneInfo = pd.DataFrame(self.expressionList.iloc[:, 0], columns=['gene_id'])

        self.sampleInfo = pd.DataFrame(range(len(self.expressionList.columns[1:])), columns=['sample_id'], index=self.expressionList.columns[1:])

        self.expressionList.index = self.expressionList.iloc[:, 0]  # gene_id
        # drop gene id columns
        self.expressionList = self.expressionList.drop([self.expressionList.columns[0]], axis=1)

    def getGeneExp(self):
        return self.expressionList

    def addGeneInfo(self, path, sep=' '):
        if not os.path.isfile(path):
            raise ValueError("file does not exist!")

        genes = pd.read_csv(path, sep=sep)
        genes.index = self.geneInfo.index
        self.geneInfo = pd.concat([self.geneInfo, genes], axis=1)

    def addSampleInfo(self, path, sep=' '):
        if not os.path.isfile(path):
            raise ValueError("file does not exist!")

        samples = pd.read_csv(path, sep=sep)
        samples.index = self.sampleInfo.index
        self.sampleInfo = pd.concat([self.sampleInfo, samples], axis=1)

    def getGeneInfo(self):
        return self.geneInfo

    def getSampleInfo(self):
        return self.sampleInfo