import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import main as WGCNA
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist


def preprocess():
    expressionList = pd.read_csv('test/input/expressionList', sep=' ')

    expressionList.index = expressionList['gene_id']
    expressionList = expressionList.drop(['gene_id'], axis=1)

    # Prepare and clean data
    # Remove rows with less than 1 TPM
    expressionList = expressionList.loc[(expressionList > 1).any(axis=1), :]

    # Check that all genes and samples have sufficiently low numbers of missing values.
    goodGenes, goodSamples, allOK = WGCNA.goodSamplesGenes(expressionList, verbose=3)
    # if not okay
    if not allOK:
        # Optionally, print the gene and sample names that were removed:
        if np.count_nonzero(goodGenes) > 0:
            print(("Removing genes:", expressionList.index[not goodGenes].values, "\n"), flush=True)
        if np.count_nonzero(goodSamples) > 0:
            print(("Removing samples:", expressionList.columns[not goodSamples].values, "\n"), flush=True)
        # Remove the offending genes and samples from the data:
        expressionList = expressionList.loc[goodGenes, goodSamples]

    print("\n\n")

    # Clustering
    sampleTree = WGCNA.hclust(pdist(expressionList.T), method="average")

    cut = 400000
    dendrogram(sampleTree, color_threshold=cut, labels=expressionList.T.index, leaf_rotation=90, leaf_font_size=8)
    plt.axhline(y=cut, c='grey', lw=1, linestyle='dashed')
    plt.title('Sample clustering to detect outliers')
    plt.xlabel('Samples')
    plt.ylabel('Distances')
    plt.tight_layout()
    plt.savefig('test/output/plots/sampleClusteringCleaning.png')

    # Determine cluster under the line
    clust = WGCNA.cutree(sampleTree, cutHeight=cut)
    # clust 0 contains the samples we want to keep.
    clust = clust.T.tolist()[0]
    index = [index for index, element in enumerate(clust) if element == 0]

    datExpr = expressionList.iloc[:, index]

    # convert trancript ID to gene ID
    for i in range(datExpr.shape[0]):
        datExpr.index.values[i] = datExpr.index[i].split(".")[0]

    datExpr = datExpr.transpose()

    datExpr.to_csv('test/output/data/data_input')


def run_WGCNA():
    datExpr = pd.read_csv('test/output/data/data_input', header=0, index_col=0)
    # Choose a set of soft-thresholding powers
    powers = list(range(1, 11)) + list(range(11, 21, 2))

    # Call the network topology analysis function
    sft = WGCNA.pickSoftThreshold(datExpr, powerVector=powers, networkType="signed", verbose=5)

    fig, ax = plt.subplots(ncols=2)
    ax[0].scatter(x=df['Gr Liv Area'], y=df['SalePrice'])
    ax[0].set_xlabel("Soft Threshold (power)")
    ax[0].set_ylabel("Scale Free Topology Model Fit,signed R^2")
    ax[0].title('Scale independence')

    ax[1].scatter(x=df['Overall Qual'], y=df['SalePrice'])
    ax[1].set_xlabel("Soft Threshold (power)")
    ax[1].set_ylabel("Mean Connectivity")
    ax[1].title('Mean connectivity')

    plt.tight_layout()
    plt.savefig('test/output/plots/summarypower.png')

    return sft


if __name__ == '__main__':
    preprocess()

    run_WGCNA()
