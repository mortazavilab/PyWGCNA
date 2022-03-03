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

    def __init__(self, geneExp=None, geneExpPath=None, sep=' '):
        if geneExpPath is not None:
            if not os.path.isfile(geneExpPath):
                raise ValueError("file does not exist!")
            else:
                self.expressionList = pd.read_csv(geneExpPath, sep=sep)
        elif geneExp is None:
            raise ValueError("geneExp and geneExpPath can not be empty at the same time!")
        else:
            if isinstance(geneExp, pd.DataFrame):
                self.expressionList = geneExp
            else:
                raise ValueError("geneExp is not data frame!")

        self.geneInfo = pd.DataFrame(self.expressionList.iloc[:, 0], columns=['gene_id'])

        self.sampleInfo = pd.DataFrame(range(len(self.expressionList.columns[1:])), columns=['sample_id'], index=self.expressionList.columns[1:])

        self.expressionList.index = self.expressionList.iloc[:, 0]  # gene_id
        # drop gene id columns
        self.expressionList = self.expressionList.drop([self.expressionList.columns[0]], axis=1)

    def addGeneInfo(self, path, sep=' '):
        if not os.path.isfile(path):
            raise ValueError("file does not exist!")

        genes = pd.read_csv(path, sep=sep)
        genes.index = self.geneInfo.index
        self.geneInfo = pd.concat([self.geneInfo, genes], axis=1)

    def updateMetadata(self, path, sep=' '):
        if not os.path.isfile(path):
            raise ValueError("file does not exist!")

        samples = pd.read_csv(path, sep=sep)
        self.sampleInfo = pd.concat([samples, self.sampleInfo],
                                    axis=1, ignore_index=True)
        self.sampleInfo = self.sampleInfo.loc[:, ~self.sampleInfo.columns.duplicated()]
