import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import main as WGCNA


def preprocess():
    expressionList = pd.read_csv('input/expressionList', header=True, index_col=0)

    # Prepare and clean data
    # Remove rows with less than 1 TPM
    expressionList = expressionList[expressionList[:, np.where(expressionList)] > 1, :]

    # Check that all genes and samples have sufficiently low numbers of missing values.
    goodGenes, goodSamples, allOK = WGCNA.goodSamplesGenes(expressionList, verbose=3)
    # if not okay
    if not allOK:
        # Optionally, print the gene and sample names that were removed:
        if np.count_nonzero(goodGenes) > 0:
            print(("Removing genes:", expressionList.index[goodGenes], "\n"), flush=True)
        if np.count_nonzero(goodSamples) > 0:
            print(("Removing samples:", expressionList.index[goodSamples], "\n"), flush=True)
        # Remove the offending genes and samples from the data:
        expressionList = expressionList[goodSamples, goodGenes]

    # Clustering
    sampleTree = WGCNA.hclust(expressionList, method="average")
    # The user should change the dimensions if the window is too large or too small.
    plt.plot(sampleTree, main="Sample clustering to detect outliers")
    # Plot a line to show the cut
    plt.axhline(y=50000, color='r', linestyle='-')
    plt.savefig('output/plots/sampleClusteringCleaning.png')

    # Determine cluster under the line
    clust = WGCNA.cutree(sampleTree, cutHeight=50000)
    # clust 1 contains the samples we want to keep.
    keepSamples = (clust == 1)

    datExpr = expressionList[keepSamples, :]
    nGenes = datExpr.shape[1]
    nSamples = datExpr.shape[2]

    # convert trancript ID to gene ID
    datExpr = datExpr.transpose()
    for i in range(datExpr.shape[1]):
        datExpr.index[i] = datExpr.index[i].split("\\.")[1]
    datExpr = datExpr.transpose()

    datExpr.to_csv('output/ata/data_input')


if __name__ == '__main__':
    print("AAA")
    preprocess()