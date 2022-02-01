import pandas as pd
import numpy as np
from scipy.stats import fisher_exact


class Comparison:
    """
            A class used to compare WGCNA to another WGCNA or any gene marker table

            Attributes
            ----------
            name1 : str
                name of first WGCNA

            name2 : str
                name of second WGCNA

            geneModule1 : dict
                gene modules of first WGCNA

            geneModule2 : dict
                gene modules of second WGCNA

            geneMarker : data frame
                gene marker of single cell data

            sc : bool
                indicate if object is WGCNA or single cell

            Methods
            -------
            compareWGCNA()
                compare two WGCNA

            compareSingleCell()
                compare WGCNA to single cell data

            """
    def __init__(self, name1="name1", name2="name2", geneModule1=None, geneModule2=None, geneMarker=None, sc=False):
        self.name1 = name1
        self.name2 = name2
        self.geneModule1 = geneModule1
        self.geneModule2 = geneModule2
        self.geneMarker = geneMarker
        self.sc = sc

        self.confusion = None

    def compareWGCNA(self):
        """
        Compare two list of modules from two bulk gene expression data set
        Returns
        -------
        compare class
        """""
        if self.name1 == self.name2:
            name1 = self.name1 + "1"
            name2 = self.name2 + "2"
        else:
            name1 = self.name1
            name2 = self.name2
        num = len(self.geneModule1.keys()) * len(self.geneModule2.keys())
        df = pd.DataFrame(columns=[name1, name2, name1 + "_size", name2 + "_size", "number", "fraction(%)", "P_value"], index=range(num))

        genes = []
        count = 0
        for i in range(len(self.geneModule1.keys())):
            node1 = self.geneModule1[self.geneModule1.keys()[i]]
            genes = genes + node1
            for j in range(len(self.geneModule2.keys())):
                node2 = self.geneModule2[self.geneModule2.keys()[j]]

                df[name1][count] = self.geneModule1.keys()[i]
                df[name2][count] = self.geneModule2.keys()[j]
                df[name1 + '_size'][count] = len(node1)
                df[name2 + '_size'][count] = len(node2)
                num = np.intersect1d(node1, node2)
                df['number'][count] = len(num)
                df['fraction(%)'][count] = len(num) / len(node2) * 100
                count = count + 1

                genes = genes + node2

        genes = list(set(genes))
        nGenes = len(genes)

        count = 0
        for i in range(len(self.geneModule1.keys())):
            for j in range(len(self.geneModule2.keys())):
                table = np.array([[nGenes - df[name1][count] - df[name2][count] + df['number'][count],
                                   df[name1][count] - df['number'][count]],
                                  [df[name2][count] - df['number'][count],
                                   df['number'][count]]])
                oddsr, p = fisher_exact(table, alternative='two-sided')
                df['P_value'][count] = p
                count = count + 1

        self.confusion = df

    def compareSingleCell(self):
        """
        Compare bulk and single cell gene expression data
        Returns
        -------
        compare class
        """""
        list_sn = np.unique(self.geneMarker['cluster'])
        num = len(self.geneModule1.keys()) * len(list_sn)
        df = pd.DataFrame(columns=["WGCNA", "sc", "WGCNA_size", "sc_size", "number", "fraction(%)", "P_value", "cellType"], index=range(num))

        genes = []
        count = 0
        for i in range(len(self.geneModule1.keys())):
            node1 = self.geneModule1[self.geneModule1.keys()[i]]
            genes = genes + node1
            for j in range(len(list_sn)):
                node2 = self.geneMarker[self.geneMarker['cluster'] == list_sn[j], :]

                df['WGCNA'][count] = self.geneModule1.keys()[i]
                df['sc'][count] = "N" + str(list_sn[j])
                df['WGCNA_size'][count] = len(node1)
                df['sc_size'][count] = len(node2)
                num = np.intersect1d(node1, node2)
                df['number'][count] = len(num)
                df['fraction(%)'][count] = len(num) / len(node2) * 100
                df['cellType'][count] = self.geneMarker['cellType'][np.where(self.geneMarker['cluster'] == list_sn[j]).tolist()[0]]
                count = count + 1

                genes = genes + node2

        genes = list(set(genes))
        nGenes = len(genes)

        count = 0
        for i in range(len(self.geneModule1.keys())):
            for j in range(len(list_sn)):
                table = np.array([[nGenes - df['WGCNA'][count] - df['sc'][count] + df['number'][count],
                                   df['WGCNA'][count] - df['number'][count]],
                                  [df['sc'][count] - df['number'][count],
                                   df['number'][count]]])
                oddsr, p = fisher_exact(table, alternative='two-sided')
                df['P_value'][count] = p
                count = count + 1

        self.confusion = df

    def getConfusionMatrix(self):
        """
        get confusion matrix as a dataframe
        Returns
        -------
        confusion matrix
        """""
        return self.confusion
