import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import pickle


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


class Comparison:
    """
    A class used to compare PyWGCNA to another PyWGCNA or any gene marker table

    :param name1: name of first WGCNA
    :type name1: str
    :param name2: name of second WGCNA
    :type name2: str
    :param geneModule1: gene modules of first WGCNA
    :type geneModule1: dict
    :param geneModule2: gene modules of second WGCNA
    :type geneModule2: dict
    :param geneMarker: gene marker of single cell data
    :type geneMarker: pandas dataframe
    :param sc: indicate if object is WGCNA or single cell
    :type sc: bool
    :param comparison: Summary of comparison results
    :type comparison: pandas dataframe

    """

    def __init__(self, name1="name1", name2="name2", geneModule1=None, geneModule2=None, geneMarker=None, sc=False):
        self.name1 = name1
        self.name2 = name2
        self.geneModule1 = geneModule1
        self.geneModule2 = geneModule2
        self.geneMarker = geneMarker
        self.sc = sc

        self.comparison = None

    def compareWGCNA(self):
        """
        Compare two list of modules from two bulk gene expression data set

        :return: update compare class that replace automatically
        :rtype: compare class
        """
        if self.name1 == self.name2:
            name1 = self.name1 + "1"
            name2 = self.name2 + "2"
        else:
            name1 = self.name1
            name2 = self.name2
        num = len(self.geneModule1.keys()) * len(self.geneModule2.keys())
        df = pd.DataFrame(columns=[name1, name2, name1 + "_size", name2 + "_size", "number", "fraction(%)", "P_value"],
                          index=range(num))

        genes = []
        count = 0
        for i in range(len(self.geneModule1.keys())):
            node1 = self.geneModule1[list(self.geneModule1.keys())[i]].gene_id.tolist()
            genes = genes + node1
            for j in range(len(self.geneModule2.keys())):
                node2 = self.geneModule2[list(self.geneModule2.keys())[j]].gene_id.tolist()

                df[name1][count] = list(self.geneModule1.keys())[i]
                df[name2][count] = list(self.geneModule2.keys())[j]
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
                table = np.array(
                    [[nGenes - df[name1 + '_size'][count] - df[name2 + '_size'][count] + df['number'][count],
                      df[name1 + '_size'][count] - df['number'][count]],
                     [df[name2 + '_size'][count] - df['number'][count],
                      df['number'][count]]])
                oddsr, p = fisher_exact(table, alternative='two-sided')
                df['P_value'][count] = p
                count = count + 1

        self.comparison = df

    def compareSingleCell(self):
        """
        Compare PyWGCNA object to single cell gene expression data

        :return: update compare class that replace automatically
        :rtype: compare class
        """
        list_sn = np.unique(self.geneMarker['cluster'])
        num = len(self.geneModule1.keys()) * len(list_sn)
        df = pd.DataFrame(
            columns=["WGCNA", "sc", "WGCNA_size", "sc_size", "number", "fraction(%)", "P_value", "cellType"],
            index=range(num))

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
                df['cellType'][count] = self.geneMarker['cellType'][
                    np.where(self.geneMarker['cluster'] == list_sn[j]).tolist()[0]]
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

        self.comparison = df

    def plotCompareWGCA(self, order1=None, order2=None, save=False):
        """
        plot comparison

        :param order1: order of modules in PyWGCNA1 you want to show in plot (name of each elements should mapped the name of modules in your first PyWGCNA)
        :type order1: list of str
        :param order2: order of modules in PyWGCNA2 you want to show in plot (name of each elements should mapped the name of modules in your second PyWGCNA)
        :type order2: list of str
        :param save: if you want to save plot as comparison.png near to your script
        :type save: bool

        """
        result = self.comparison.copy(deep=True)
        result['-log10(P_value)'] = -1 * np.log10(result['P_value'].astype(np.float64))

        if self.name1 == self.name2:
            name1 = self.name1 + "1"
            name2 = self.name2 + "2"
        else:
            name1 = self.name1
            name2 = self.name2

        result.drop(labels=np.where(result[name1] == 'grey')[0].tolist(),
                    axis=0,
                    inplace=True)
        result.reset_index(drop=True, inplace=True)
        result.drop(labels=np.where(result[name2] == 'grey')[0].tolist(),
                    axis=0,
                    inplace=True)
        result.reset_index(drop=True, inplace=True)

        result.loc[np.where(result['fraction(%)'] == 0)[0].tolist(), 'fraction(%)'] = np.nan
        result.loc[np.where(result['fraction(%)'] == 0)[0].tolist(), 'fraction(%)'] = np.nan

        result.loc[np.isinf(result['-log10(P_value)']), '-log10(P_value)'] = np.max(
            result['-log10(P_value)'][np.isfinite(result['-log10(P_value)'])]) + 1

        grey = result.copy(deep=True)
        result.loc[np.where(result['P_value'] > 0.01)[0].tolist(), '-log10(P_value)'] = np.nan

        result.dropna(axis=0, inplace=True)
        result.reset_index(drop=True, inplace=True)

        grey.loc[np.where(grey['P_value'] <= 0.01)[0].tolist(), '-log10(P_value)'] = np.nan
        grey.dropna(axis=0, inplace=True)
        grey.reset_index(drop=True, inplace=True)

        if order1 is not None:
            result[name1] = pd.Categorical(result[name1], order1)
            result.sort_values(by=[name1], inplace=True)

            grey[name1] = pd.Categorical(grey[name1], order1)
            grey.sort_values(by=[name1], inplace=True)

        if order2 is not None:
            result[name2] = pd.Categorical(result[name2], order2)
            result.sort_values(by=[name2], inplace=True)

            grey[name2] = pd.Categorical(grey[name2], order2)
            grey.sort_values(by=[name2], inplace=True)

        fig, ax = plt.subplots(figsize=(max(7, len(np.unique(result[name1])) / 3),
                                        max(5, len(np.unique(result[name2])) / 3)))
        scatter = ax.scatter(x=result[name1],
                             y=result[name2],
                             s=result['fraction(%)'].astype(float)*4,
                             c=result['-log10(P_value)'],
                             alpha=0.8,
                             cmap='viridis_r',
                             vmin=np.min(result['fraction(%)']),
                             vmax=np.max(result['fraction(%)']))
        # Add a colorbar
        fig.colorbar(scatter, shrink=0.25, label='-log10(P_value)')

        ax.scatter(x=grey[name1],
                   y=grey[name2],
                   s=grey['fraction(%)'].astype(float)*4,
                   c='grey',
                   alpha=0.8,
                   vmin=np.min(grey['fraction(%)']),
                   vmax=np.max(grey['fraction(%)']))

        # produce a legend with the unique colors from the scatter
        kw = dict(prop="sizes", num=4, color='black', fmt="{x:.1f} %",
                  func=lambda s: s/4)
        legend2 = ax.legend(*scatter.legend_elements(**kw),
                            bbox_to_anchor=(1.05, 0.98), loc="upper left",
                            title="fraction(%)")

        plt.xticks(rotation=90)
        plt.xlabel(name1 + " modules")
        plt.ylabel(name2 + " modules")

        plt.tight_layout()

        if save:
            plt.savefig('comparison_' + name1 + '_' + name2 + '.png')
        plt.show()

    def saveComparison(self):
        """
        save comparison object as comparison.p near to the script
        
        """
        print(f"{BOLD}{OKBLUE}Saving comparison as comparison.p{ENDC}")

        picklefile = open('comparison.p', 'wb')
        pickle.dump(self, picklefile)
        picklefile.close()

