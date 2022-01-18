import math
import numpy as np
import pandas as pd
import scipy.stats as stats
import statistics
import sys
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, cut_tree, dendrogram
import networkx as nx
from statsmodels.formula.api import ols
import resource
from matplotlib import colors as mcolors
from sklearn.impute import KNNImputer
from sklearn.preprocessing import scale
import random

random.seed(10)

# remove runtime warning (divided by zero)
np.seterr(divide='ignore', invalid='ignore')


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# public values
networkTypes = ["unsigned", "signed", "signed hybrid"]
adjacencyTypes = ["unsigned", "signed", "signed hybrid"]
TOMTypes = ["NA", "unsigned", "signed"]
TOMDenoms = ["min", "mean"]


def replaceMissing(x, replaceWith):
    if replaceWith:
        if x.isnumeric():
            replaceWith = 0
        elif x.isalpha():
            replaceWith = ""
        else:
            sys.exit("Need 'replaceWith'.")

    x = x.fillna(replaceWith)
    return x


def checkAndScaleWeights(weights, expr, scaleByMax=True, verbose=1):
    if weights is None:
        return weights

    weights = np.asmatrix(weights)
    if expr.shape != weights.shape:
        sys.exit("When 'weights' are given, they must have the same dimensions as 'expr'.")
    if (weights < 0).any():
        sys.exit("Found negative weights. All weights must be non-negative.")

    nf = np.isinf(weights)
    if any(nf):
        if verbose > 0:
            print(
                f"{bcolors.WARNING}Found non-finite weights. The corresponding data points will be removed.{bcolors.ENDC}")
            weights[nf] = None

    if scaleByMax:
        maxw = np.amax(weights, axis=0)
        maxw[maxw == 0] = 1
        weights = weights / np.reshape(maxw, weights.shape)

    return weights


# Filter genes with too many missing entries
def goodGenesFun(datExpr, weights=None, useSamples=None, useGenes=None, minFraction=1 / 2,
                 minNSamples=4, minNGenes=4, tol=None, minRelativeWeight=0.1, verbose=1):
    if not datExpr.apply(lambda s: pd.to_numeric(s, errors='coerce').notnull().all()).all():
        sys.exit("datExpr must contain numeric data.")

    weights = checkAndScaleWeights(weights, datExpr, scaleByMax=True)

    if tol is None:
        tol = 1e-10 * datExpr.abs().max().max()
    if useGenes is None:
        useGenes = np.repeat(True, datExpr.shape[0])
    if useSamples is None:
        useSamples = np.repeat(True, datExpr.shape[1])

    if len(useGenes) != datExpr.shape[0]:
        sys.exit("Length of nGenes is not compatible with number of columns in datExpr.")
    if len(useSamples) != datExpr.shape[1]:
        sys.exit("Length of nSamples is not compatible with number of rows in datExpr.")

    nSamples = sum(useSamples)
    nGenes = sum(useGenes)
    if weights is None:
        nPresent = datExpr.loc[useGenes, useSamples].notna().sum(axis=1)
    else:
        nPresent = (datExpr.loc[useGenes, useSamples].notna() and
                    weights.loc[useGenes, useSamples] > minRelativeWeight).sum(axis=1)

    gg = useGenes
    gg[np.logical_and(useGenes, nPresent < minNSamples)] = False

    if weights is None:
        var = np.var(datExpr.loc[gg, useSamples], axis=1)
        # var = var.sort_index(inplace=True)
    else:
        # need to be fix
        # TODO:colWeightedVars
        var = np.var(datExpr, w=weights)

    var[np.isnan(var)] = 0
    nNAsGenes = datExpr.loc[gg, useSamples].isna().sum(axis=1)
    gg[gg] = np.logical_and(np.logical_and(nNAsGenes < (1 - minFraction) * nSamples, var > tol ** 2),
                            nSamples - nNAsGenes >= minNSamples)

    if sum(gg) < minNGenes:
        sys.exit("Too few genes with valid expression levels in the required number of samples.")
    if verbose > 0 and nGenes - sum(gg) > 0:
        print("\n\n  ..Excluding", nGenes - sum(gg),
              "genes from the calculation due to too many missing samples or zero variance.\n\n", flush=True)

    return gg


# Filter samples with too many missing entries
def goodSamplesFun(datExpr, weights=None, useSamples=None, useGenes=None, minFraction=1 / 2,
                   minNSamples=4, minNGenes=4, minRelativeWeight=0.1, verbose=1):
    if useGenes is None:
        useGenes = np.repeat(True, datExpr.shape[0])

    if useSamples is None:
        useSamples = np.repeat(True, datExpr.shape[1])

    if len(useGenes) != datExpr.shape[0]:
        sys.exit("Length of nGenes is not compatible with number of columns in datExpr.")

    if len(useSamples) != datExpr.shape[1]:
        sys.exit("Length of nSamples is not compatible with number of rows in datExpr.")

    weights = checkAndScaleWeights(weights, datExpr, scaleByMax=True)
    nSamples = sum(useSamples)
    nGenes = sum(useGenes)
    if weights is None:
        nNAsSamples = np.sum((datExpr.loc[useGenes, useSamples]).isnull(), axis=0)
    else:
        nNAsSamples = np.sum(np.logical_or(datExpr[useGenes, useSamples],
                                           replaceMissing(weights[useGenes, useSamples] < minRelativeWeight, True))
                             .isnull(), axis=0)

    goodSamples = useSamples
    goodSamples[useSamples] = np.logical_and((nNAsSamples < (1 - minFraction) * nGenes),
                                             (nGenes - nNAsSamples >= minNGenes))

    if sum(goodSamples) < minNSamples:
        sys.exit("Too few samples with valid expression levels for the required number of genes.")

    if verbose > 0 and (nSamples - sum(goodSamples) > 0):
        print("  ..Excluding", nSamples - sum(goodSamples),
              "samples from the calculation due to too many missing genes.", flush=True)

    return goodSamples


# Check that all genes and samples have sufficiently low numbers of missing values.
def goodSamplesGenes(datExpr, weights=None, minFraction=1 / 2, minNSamples=4, minNGenes=4, tol=None,
                     minRelativeWeight=0.1, verbose=1):
    goodGenes = None
    goodSamples = None
    nBadGenes = 0
    nBadSamples = 0
    changed = True
    iter = 1
    if verbose > 0:
        print("Flagging genes and samples with too many missing values...", flush=True)
    while changed:
        if verbose > 0:
            print(" ..step", iter, flush=True)
        goodGenes = goodGenesFun(datExpr, weights, goodSamples, goodGenes, minFraction=minFraction,
                                 minNSamples=minNSamples, minNGenes=minNGenes, minRelativeWeight=minRelativeWeight,
                                 tol=tol, verbose=verbose - 1)
        goodSamples = goodSamplesFun(datExpr, weights, goodSamples, goodGenes, minFraction=minFraction,
                                     minNSamples=minNSamples, minNGenes=minNGenes, minRelativeWeight=minRelativeWeight,
                                     verbose=verbose - 1)
        changed = np.logical_or((np.logical_not(goodGenes).sum() > nBadGenes),
                                (np.logical_not(goodSamples).sum() > nBadSamples))
        nBadGenes = np.logical_not(goodGenes).sum()
        nBadSamples = np.logical_not(goodSamples).sum()
        iter = iter + 1

    allOK = (nBadGenes + nBadSamples == 0)

    return goodGenes, goodSamples, allOK


def hclust(d, method="complete"):
    METHODS = ["single", "complete", "average", "weighted", "centroid"]

    if method not in METHODS:
        sys.exit("Invalid clustering method.")

    if method == -1:
        sys.exit("Ambiguous clustering method.")

    dendrogram = linkage(d, method=method)

    return dendrogram


# Determine cluster under the line
def cutree(sampleTree, cutHeight=50000.0):
    cutTree = cut_tree(sampleTree, height=cutHeight)

    return cutTree


def checkSimilarity(adjMat, min=0, max=1):
    dim = adjMat.shape
    if dim is None or len(dim) != 2:
        sys.exit("adjacency is not two-dimensional")

    if not (all(ele.isdigit() for ele in adjMat)):
        sys.exit("adjacency is not numeric")

    if dim[1] != dim[2]:
        sys.exit("adjacency is not square")

    if np.max(np.abs(adjMat - adjMat.transpose())) > 1e-12:
        sys.exit("adjacency is not symmetric")

    if np.min(adjMat) < min or np.max(adjMat) > max:
        sys.exit(("some entries are not between", min, "and", max))


def calBlockSize(matrixSize, rectangularBlocks=True, maxMemoryAllocation=None, overheadFactor=3):
    if maxMemoryAllocation is None:
        maxAlloc = resource.getrlimit(resource.RLIMIT_AS)[1]
    else:
        maxAlloc = maxMemoryAllocation / 8

    maxAlloc = maxAlloc / overheadFactor

    if rectangularBlocks:
        blockSz = math.floor(maxAlloc / matrixSize)
    else:
        blockSz = math.floor(math.sqrt(maxAlloc))

    return min(matrixSize, blockSz)


# Calculation of fitting statistics for evaluating scale free topology fit.
def scaleFreeFitIndex(k, nBreaks=10):
    df = pd.DataFrame({'data': k})
    df['discretized_k'] = pd.cut(df['data'], nBreaks)
    dk = df.groupby('discretized_k').mean()  # tapply(k, discretized_k, mean)
    dk = pd.DataFrame(dk.reset_index())
    dk.columns = ['discretized_k', 'dk']
    p_dk = df['discretized_k'].value_counts() / len(k)  # as.vector(tapply(k, discretized.k, length)/length(k))
    p_dk = pd.DataFrame(p_dk.reset_index())
    p_dk.columns = ['discretized_k', 'p_dk']
    breaks1 = np.linspace(start=min(k), stop=max(k), num=nBreaks + 1)
    y, edges = np.histogram(df['data'], bins=breaks1)
    dk2 = 0.5 * (edges[1:] + edges[:-1])
    df = pd.merge(dk, p_dk, on='discretized_k')
    if df['dk'].isnull().values.any():
        df.loc[df['dk'].isnull().values, 'dk'] = dk2[df['dk'].isnull().values]
    if np.any(df['dk'] == 0):
        df.loc[df['dk'] == 0, 'dk'] = dk2[df['dk'] == 0]
    if df['p_dk'].isnull().values.any():
        df.loc[df['p_dk'].isnull().values, 'p_dk'] = 0
    df['log_dk'] = np.log10(df['dk'])
    df['log_p_dk'] = np.log10(df['p_dk'] + 1e-09)
    df['log_p_dk_10'] = np.power(10, df['log_dk'])

    model1 = ols(formula='log_p_dk ~ log_dk', data=df).fit()
    model2 = ols(formula='log_p_dk ~ log_dk + log_p_dk_10', data=df).fit()
    dfout = pd.DataFrame({'Rsquared.SFT': [model1.rsquared],
                          'slope.SFT': [model1.params.values[1]],
                          'truncatedExponentialAdjRsquared': [model2.rsquared_adj]})
    return dfout


# Call the network topology analysis function
def pickSoftThreshold(data, dataIsExpr=True, weights=None, RsquaredCut=0.9,
                      powerVector=np.concatenate((np.arange(1, 10, 1), np.arange(12, 20, 2))),
                      nBreaks=10, blockSize=None, corOptions=None, networkType="unsigned",
                      moreNetworkConcepts=False, gcInterval=None):
    powerVector = np.sort(powerVector)
    intType = networkTypes.index(networkType)
    if intType is None:
        sys.exit(("Unrecognized 'networkType'. Recognized values are", str(networkTypes)))

    nGenes = data.shape[1]
    if nGenes < 3:
        sys.exit("The input data data contain fewer than 3 rows (nodes).\n"
                 "This would result in a trivial correlation network.")

    if not dataIsExpr:
        checkSimilarity(data)
        if any(np.diag(data) != 1):
            data = np.where(np.diag(data), 1)

    if blockSize is None:
        blockSize = calBlockSize(nGenes, rectangularBlocks=True, maxMemoryAllocation=2 ** 30)
        print("pickSoftThreshold: will use block size ", blockSize, flush=True)

    if gcInterval is None or len(gcInterval) == 0:
        gcInterval = 4 * blockSize

    colname1 = ["Power", "SFT.R.sq", "slope", "truncated R.sq", "mean(k)", "median(k)", "max(k)"]

    if moreNetworkConcepts:
        colname1 = colname1.append(["Density", "Centralization", "Heterogeneity"])

    datout = pd.DataFrame(np.full((len(powerVector), len(colname1)), 666), columns=colname1, dtype=object)
    datout['Power'] = powerVector
    print("pickSoftThreshold: calculating connectivity for given powers...", flush=True)

    datk = np.zeros((nGenes, len(powerVector)))
    nPowers = len(powerVector)
    startG = 0
    lastGC = 0
    if corOptions is None:
        corOptions = pd.DataFrame()
        corOptions['x'] = [data]
    else:
        corOptions['x'] = [data]

    if weights is not None:
        if not dataIsExpr:
            sys.exit("Weights can only be used when 'data' represents expression data ('dataIsExpr' must be TRUE).")
        if data.shape != weights.shape:
            sys.exit("When 'weights' are given, dimensions of 'data' and 'weights' must be the same.")
        corOptions['weights.x'] = weights

    while startG < nGenes:
        endG = min(startG + blockSize, nGenes)
        print("\t..working on genes", (startG + 1), "through", endG, "of", nGenes, flush=True)

        useGenes = list(range(startG, endG))
        nGenes1 = len(useGenes)
        if dataIsExpr:
            corOptions['y'] = [data.iloc[:, useGenes]]
            if weights is not None:
                corOptions['weights.y'] = [weights.iloc[:, useGenes]]
            corx = np.corrcoef(corOptions.x[0], corOptions.y[0], rowvar=False)
            corx = corx[0:corOptions.x[0].shape[1], useGenes]
            if intType == 0:
                corx = abs(corx)
            elif intType == 1:
                corx = (1 + corx) / 2
            elif intType == 2:
                corx[corx < 0] = 0

            if np.count_nonzero(np.isnan(corx)) != 0:
                print(f"{bcolors.WARNING}Some correlations are NA in block {str(startG)} : {str(endG)}.{bcolors.ENDC}")
        else:
            corx = data[:, useGenes]

        corx[useGenes, list(range(len(useGenes)))] = 1
        datk_local = np.empty((nGenes1, nPowers))
        datk_local[:] = np.nan
        corxPrev = np.ones(corx.shape)
        powerVector1 = [0]
        powerVector1.extend(powerVector[:-1])
        powerSteps = powerVector - powerVector1
        uniquePowerSteps = np.unique(powerSteps)

        def func(power):
            return corx ** power

        corxPowers = pd.DataFrame()
        for p in uniquePowerSteps:
            corxPowers[p] = [func(p)]

        for j in range(nPowers):
            corxCur = corxPrev * corxPowers[powerSteps[j]][0]
            datk_local[:, j] = np.nansum(corxCur, axis=0) - 1
            corxPrev = corxCur

        datk[startG:endG, :] = datk_local

        startG = endG
        if 0 < gcInterval < startG - lastGC:
            lastGC = startG

    pd.DataFrame(datk).to_csv('test/output/data/datk')
    print(datk.shape)

    for i in range(len(powerVector)):
        khelp = datk[:, i]
        SFT1 = scaleFreeFitIndex(k=khelp, nBreaks=nBreaks)
        datout.loc[i, 'SFT.R.sq'] = SFT1.loc[0, 'Rsquared.SFT']
        datout.loc[i, 'slope'] = SFT1.loc[0, 'slope.SFT']
        datout.loc[i, 'truncated R.sq'] = SFT1.loc[0, 'truncatedExponentialAdjRsquared']
        datout.loc[i, 'mean(k)'] = statistics.mean(khelp)
        datout.loc[i, 'median(k)'] = statistics.median(khelp)
        datout.loc[i, 'max(k)'] = max(khelp)

        if moreNetworkConcepts:
            Density = sum(khelp) / (nGenes * (nGenes - 1))
            datout.loc[i, 'Density'] = Density
            Centralization = nGenes * (max(khelp) - statistics.mean(khelp)) / ((nGenes - 1) * (nGenes - 2))
            datout.loc[i, 'Centralization'] = Centralization
            Heterogeneity = np.sqrt(nGenes * sum(khelp ^ 2) / sum(khelp) ^ 2 - 1)
            datout.loc[i, 'Heterogeneity'] = Heterogeneity

    print(datout)

    # detect threshold more than 0.9 by default
    ind = datout['SFT.R.sq'] > RsquaredCut
    if np.sum(ind) > 0:
        powerEstimate = np.min(powerVector[ind])
        print("Selected power to have scale free network is ", powerEstimate, "\n\n")
    else:
        ind = np.argsort(datout['SFT.R.sq']).tolist()
        powerEstimate = powerVector[ind[-1]]
        print("No power detected to have scale free network!\n",
              "Found the best given power which is ", powerEstimate, "\n\n")

    return powerEstimate, datout


def adjacency(datExpr, selectCols=None, adjacencyType="unsigned", power=6, corOptions=pd.DataFrame(), weights=None,
              weightArgNames=None):
    if weightArgNames is None:
        weightArgNames = ["weights.x", "weights.y"]
    intType = adjacencyTypes.index(adjacencyType)
    if intType is None:
        sys.exit(("Unrecognized 'type'. Recognized values are", str(adjacencyTypes)))
    weights = checkAndScaleWeights(weights, datExpr, scaleByMax=False)
    if weights is not None:
        if selectCols is None:
            if isinstance(corOptions, pd.DataFrame):
                weightOpt = pd.DataFrame({'weights.x': weights})
                weightOpt.index = weightArgNames[0]
            else:
                weightOpt = weightArgNames[0] + " = weights"
        else:
            if isinstance(corOptions, pd.DataFrame):
                weightOpt = pd.DataFrame({'weights.x': weights, 'weights.y': weights[:, selectCols]})
                weightOpt.index = weightArgNames[0:2]
            else:
                weightOpt = weightArgNames[1] + " = weights, " + weightArgNames[2] + " = weights[, selectCols]"
    else:
        if isinstance(corOptions, pd.DataFrame):
            weightOpt = pd.DataFrame()
        else:
            weightOpt = ""

    if selectCols is None:
        cor_mat = np.corrcoef(datExpr.T)  # cor_mat = do.call(corFnc.fnc, c(list(x = datExpr), weightOpt, corOptions))
    else:
        cor_mat = np.corrcoef(x=datExpr, y=datExpr[:, selectCols])  # , weightOpt, corOptions)

    if intType == 0:
        cor_mat = abs(cor_mat)
    elif intType == 1:
        cor_mat = (1 + cor_mat) / 2
    elif intType == 2:
        cor_mat[cor_mat < 0] = 0

    return cor_mat ** power


def checkAdjMat(adjMat, min=0, max=1):
    shape = adjMat.shape
    if shape is None or len(shape) != 2:
        sys.exit("adjacency is not two-dimensional")

    if not issubclass(adjMat.dtype.type, np.floating):
        sys.exit("adjacency is not numeric")
    if shape[0] != shape[1]:
        sys.exit("adjacency is not square")
    if np.max(np.fabs(np.subtract(adjMat, adjMat.T))) > 1e-12:
        sys.exit("adjacency is not symmetric")
    if np.min(adjMat) < min or np.max(adjMat) > max:
        sys.exit(("some entries are not between", min, "and", max))


def TomSimilarityFromAdj(adjMat, TOMDenom, TOMType):
    # Prepare adjacency
    np.fill_diagonal(adjMat, 0)
    # Prepare TOM
    tom = np.zeros_like(adjMat)
    # Compute TOM
    L = np.matmul(adjMat, adjMat)
    ki = adjMat.sum(axis=1)
    kj = adjMat.sum(axis=0)
    if TOMDenom == 0:  # min
        MINK = np.array([np.minimum(ki_, kj) for ki_ in ki])
    else:  # mean
        MINK = np.array([((ki_ + kj) / 2) for ki_ in ki])
    if TOMType == 0:  # unsigned
        tom = (L + adjMat) / (MINK + 1 - adjMat)
    else:  # signed
        tom = np.fabs((L + adjMat)) / (MINK + 1 - np.fabs(adjMat))
    np.fill_diagonal(tom, 1)
    return tom


def TOMsimilarity(adjMat, TOMType="unsigned", TOMDenom="min"):
    TOMTypeC = TOMTypes.index(TOMType)
    if TOMTypeC is None:
        sys.exit(("Invalid 'TOMType'. Recognized values are", str(TOMTypes)))
    if TOMTypeC == 0:
        sys.exit("'TOMType' cannot be 'none' for this function.")
    TOMDenomC = TOMDenoms.index(TOMDenom)
    if TOMDenomC is None:
        sys.exit(("Invalid 'TOMDenom'. Recognized values are", str(TOMDenoms)))

    min = 0
    if TOMTypeC == 2:
        min = -1
    checkAdjMat(adjMat, min=min, max=1)
    np.nan_to_num(adjMat, copy=False, nan=0)

    print("..connectivity..")

    tom = TomSimilarityFromAdj(adjMat, TOMDenomC, TOMTypeC)

    print("..done..\n\n")

    return pd.dataframe(tom)


def interpolate(data, index):
    i = round(index)
    n = len(data)
    if i < 1:
        return data[1]
    if i >= n:
        return data[n]
    r = index - i
    return data[i] * (1 - r) + data[i + 1] * r


def coreSizeFunc(BranchSize, minClusterSize):
    BaseCoreSize = minClusterSize / 2 + 1
    if BaseCoreSize < BranchSize:
        CoreSize = int(BaseCoreSize + math.sqrt(BranchSize - BaseCoreSize))
    else:
        CoreSize = BranchSize

    return CoreSize


def cutreeHybrid(dendro, distM, cutHeight=None, minClusterSize=20, deepSplit=1,
                 maxCoreScatter=None, minGap=None, maxAbsCoreScatter=None,
                 minAbsGap=None, minSplitHeight=None, minAbsSplitHeight=None,
                 externalBranchSplitFnc=None, nExternalSplits=0, minExternalSplit=None,
                 externalSplitOptions=pd.DataFrame(), externalSplitFncNeedsDistance=None,
                 assumeSimpleExternalSpecification=True, pamStage=True,
                 pamRespectsDendro=True, useMedoids=False, maxPamDist=None,
                 respectSmallClusters=True, verbose=2):
    tmp = dendro[:, 0] > dendro.shape[0]
    dendro[tmp, 0] = dendro[tmp, 0] - dendro.shape[0]
    dendro[np.logical_not(tmp), 0] = -1 * (dendro[np.logical_not(tmp), 0] + 1)
    tmp = dendro[:, 1] > dendro.shape[0]
    dendro[tmp, 1] = dendro[tmp, 1] - dendro.shape[0]
    dendro[np.logical_not(tmp), 1] = -1 * (dendro[np.logical_not(tmp), 1] + 1)

    chunkSize = dendro.shape[0]

    if maxPamDist is None:
        maxPamDist = cutHeight

    nMerge = dendro.shape[0]
    if nMerge < 1:
        sys.exit("The given dendrogram is suspicious: number of merges is zero.")
    if distM is None:
        sys.exit("distM must be non-NULL")
    if distM.shape is None:
        sys.exit("distM must be a matrix.")
    if distM.shape[0] != nMerge + 1 or distM.shape[1] != nMerge + 1:
        sys.exit("distM has incorrect dimensions.")
    if pamRespectsDendro and not respectSmallClusters:
        print("cutreeHybrid Warning: parameters pamRespectsDendro (TRUE) "
              "and respectSmallClusters (FALSE) imply contradictory intent.\n"
              "Although the code will work, please check you really intented "
              "these settings for the two arguments.", flush=True)
    if any(np.diag(distM) != 0):
        np.fill_diagonal(distM, 0)
    refQuantile = 0.05
    refMerge = round(nMerge * refQuantile) - 1
    if refMerge < 0:
        refMerge = 0
    refHeight = dendro[refMerge, 2]
    if cutHeight is None:
        cutHeight = 0.99 * (np.max(dendro[:, 2]) - refHeight) + refHeight
        if verbose > 0:
            print("..cutHeight not given, setting it to", round(cutHeight, 3),
                  " ===>  99% of the (truncated) height range in dendro.", flush=True)
    else:
        if cutHeight > np.max(dendro[:, 2]):
            cutHeight = np.max(dendro[:, 2])
    if maxPamDist is None:
        maxPamDist = cutHeight
    nMergeBelowCut = np.count_nonzero(dendro[:, 2] <= cutHeight)
    if nMergeBelowCut < minClusterSize:
        if verbose > 0:
            print("cutHeight set too low: no merges below the cut.", flush=True)
        return pd.DataFrame({'labels': np.repeat(0, nMerge + 1, axis=0)})

    if externalBranchSplitFnc is not None:
        nExternalSplits = len(externalBranchSplitFnc)
        if len(minExternalSplit) < 1:
            sys.exit("'minExternalBranchSplit' must be given.")
        if assumeSimpleExternalSpecification and nExternalSplits == 1:
            externalSplitOptions = pd.DataFrame(externalSplitOptions)
        # TODO: externalBranchSplitFnc = lapply(externalBranchSplitFnc, match.fun)
        for es in range(nExternalSplits):
            externalSplitOptions['tree'][es] = dendro
            if len(externalSplitFncNeedsDistance) == 0 or externalSplitFncNeedsDistance[es]:
                externalSplitOptions['dissimMat'][es] = distM

    MxBranches = nMergeBelowCut
    branch_isBasic = np.repeat(True, MxBranches, axis=0)
    branch_isTopBasic = np.repeat(True, MxBranches, axis=0)
    branch_failSize = np.repeat(False, MxBranches, axis=0)
    branch_rootHeight = np.repeat(np.nan, MxBranches, axis=0)
    branch_size = np.repeat(2, MxBranches, axis=0)
    branch_nMerge = np.repeat(1, MxBranches, axis=0)
    branch_nSingletons = np.repeat(2, MxBranches, axis=0)
    branch_nBasicClusters = np.repeat(0, MxBranches, axis=0)
    branch_mergedInto = np.repeat(0, MxBranches, axis=0)
    branch_attachHeight = np.repeat(np.nan, MxBranches, axis=0)
    branch_singletons = pd.DataFrame()
    branch_basicClusters = pd.DataFrame()
    branch_mergingHeights = pd.DataFrame()
    branch_singletonHeights = pd.DataFrame()
    nBranches = -1

    defMCS = [0.64, 0.73, 0.82, 0.91, 0.95]
    defMG = [(1.0 - defMC) * 3.0 / 4.0 for defMC in defMCS]
    nSplitDefaults = len(defMCS)
    if isinstance(deepSplit, bool):
        deepSplit = pd.to_numeric(deepSplit) * (nSplitDefaults - 2)
    if deepSplit < 0 or deepSplit > nSplitDefaults:
        msg = "Parameter deepSplit (value" + str(deepSplit) + \
              ") out of range: allowable range is 0 through", str(nSplitDefaults - 1)
        sys.exit(msg)
    if maxCoreScatter is None:
        maxCoreScatter = interpolate(defMCS, deepSplit)
    if minGap is None:
        minGap = interpolate(defMG, deepSplit)
    if maxAbsCoreScatter is None:
        maxAbsCoreScatter = refHeight + maxCoreScatter * (cutHeight - refHeight)
    if minAbsGap is None:
        minAbsGap = minGap * (cutHeight - refHeight)
    if minSplitHeight is None:
        minSplitHeight = 0
    if minAbsSplitHeight is None:
        minAbsSplitHeight = refHeight + minSplitHeight * (cutHeight - refHeight)
    nPoints = nMerge + 1
    IndMergeToBranch = np.repeat(-1, nMerge, axis=0)
    onBranch = np.repeat(0, nPoints, axis=0)
    RootBranch = 0
    if verbose > 2:
        print("..Going through the merge tree", flush=True)

    mergeDiagnostics = pd.DataFrame({'smI': np.repeat(np.nan, nMerge, axis=0),
                                     'smSize': np.repeat(np.nan, nMerge, axis=0),
                                     'smCrSc': np.repeat(np.nan, nMerge, axis=0),
                                     'smGap': np.repeat(np.nan, nMerge, axis=0),
                                     'lgI': np.repeat(np.nan, nMerge, axis=0),
                                     'lgSize': np.repeat(np.nan, nMerge, axis=0),
                                     'lgCrSc': np.repeat(np.nan, nMerge, axis=0),
                                     'lgGap': np.repeat(np.nan, nMerge, axis=0),
                                     'merged': np.repeat(np.nan, nMerge, axis=0)})
    if externalBranchSplitFnc is not None:
        externalMergeDiags = pd.DataFrame(np.nan, index=list(range(nMerge)), columns=list(range(nExternalSplits)))

    extender = np.repeat(0, chunkSize, axis=0)

    for merge in range(nMerge):
        if dendro[merge, 2] <= cutHeight:
            if dendro[merge, 0] < 0 and dendro[merge, 1] < 0:
                nBranches = nBranches + 1
                branch_isBasic[nBranches] = True
                branch_isTopBasic[nBranches] = True
                branch_singletons.insert(nBranches, nBranches,
                                         np.concatenate((-1 * dendro[merge, 0:2], extender), axis=0))
                branch_basicClusters.insert(nBranches, nBranches, extender)
                branch_mergingHeights.insert(nBranches, nBranches,
                                             np.concatenate((np.repeat(dendro[merge, 2], 2), extender), axis=0))
                branch_singletonHeights.insert(nBranches, nBranches,
                                               np.concatenate((np.repeat(dendro[merge, 2], 2), extender), axis=0))
                IndMergeToBranch[merge] = nBranches
                RootBranch = nBranches
            elif np.sign(dendro[merge, 0]) * np.sign(dendro[merge, 1]) < 0:
                clust = IndMergeToBranch[int(np.max(dendro[merge, 0:2])) - 1]

                if clust == -1:
                    sys.exit("Internal error: a previous merge has no associated cluster. Sorry!")

                gene = -1 * int(np.min(dendro[merge, 0:2]))
                ns = branch_nSingletons[clust]
                nm = branch_nMerge[clust]

                if branch_isBasic[clust]:
                    branch_singletons.loc[ns, clust] = gene
                    branch_singletonHeights.loc[ns, clust] = dendro[merge, 2]
                else:
                    onBranch[int(gene)] = clust

                branch_mergingHeights.loc[nm, clust] = dendro[merge, 2]
                branch_size[clust] = branch_size[clust] + 1
                branch_nMerge[clust] = nm + 1
                branch_nSingletons[clust] = ns + 1
                IndMergeToBranch[merge] = clust
                RootBranch = clust
            else:
                clusts = IndMergeToBranch[dendro[merge, 0:2].astype(int) - 1]
                sizes = branch_size[clusts]
                rnk = np.argsort(sizes)
                small = clusts[rnk[0]]
                large = clusts[rnk[1]]
                sizes = sizes[rnk]

                if branch_isBasic[small]:
                    coresize = coreSizeFunc(branch_nSingletons[small], minClusterSize) - 1
                    Core = branch_singletons.loc[0:coresize, small] - 1
                    Core = Core.astype(int).tolist()
                    SmAveDist = np.mean(distM.iloc[Core, Core].sum() / coresize)
                else:
                    SmAveDist = 0

                if branch_isBasic[large]:
                    coresize = coreSizeFunc(branch_nSingletons[large], minClusterSize) - 1
                    Core = branch_singletons.loc[0:coresize, large] - 1
                    Core = Core.astype(int).tolist()
                    LgAveDist = np.mean(distM.iloc[Core, Core].sum() / coresize)
                else:
                    LgAveDist = 0

                mergeDiagnostics.loc[merge, :] = [small, branch_size[small], SmAveDist, dendro[merge, 2] - SmAveDist,
                                                  large, branch_size[large], LgAveDist, dendro[merge, 2] - LgAveDist,
                                                  None]
                SmallerScores = [branch_isBasic[small], branch_size[small] < minClusterSize,
                                 SmAveDist > maxAbsCoreScatter, dendro[merge, 2] - SmAveDist < minAbsGap,
                                 dendro[merge, 2] < minAbsSplitHeight]
                if SmallerScores[0] * np.count_nonzero(SmallerScores[1:]) > 0:
                    DoMerge = True
                    SmallerFailSize = not (SmallerScores[2] | SmallerScores[3])
                else:
                    LargerScores = [branch_isBasic[large],
                                    branch_size[large] < minClusterSize, LgAveDist > maxAbsCoreScatter,
                                    dendro[merge, 2] - LgAveDist < minAbsGap,
                                    dendro[merge, 2] < minAbsSplitHeight]
                    if LargerScores[0] * np.count_nonzero(LargerScores[1:]) > 0:
                        DoMerge = True
                        SmallerFailSize = not (LargerScores[2] | LargerScores[3])
                        x = small
                        small = large
                        large = x
                        sizes = sizes.reverse()
                    else:
                        DoMerge = False

                if DoMerge:
                    mergeDiagnostics['merged'][merge] = 1

                if not DoMerge and nExternalSplits > 0 and branch_isBasic[small] and branch_isBasic[large]:
                    if verbose > 4:
                        print("Entering external split code on merge ", merge, flush=True)
                    branch1 = branch_singletons[[large]][0:sizes[1]]
                    branch2 = branch_singletons[[small]][0:sizes[0]]
                    if verbose > 4:
                        print("  ..branch lengths: ", sizes[0], ", ", sizes[1], flush=True)
                    es = 0
                    while es < nExternalSplits and not DoMerge:
                        es = es + 1
                        args = pd.DataFrame({'externalSplitOptions': externalSplitOptions[[es]],
                                             'branch1': branch1, 'branch2': branch2})
                        # TODO: extSplit = do.call(externalBranchSplitFnc[[es]], args)
                        extSplit = None
                        DoMerge = extSplit < minExternalSplit[es]
                        externalMergeDiags[merge, es] = extSplit
                        mergeDiagnostics['merged'][merge] = 0
                        if DoMerge:
                            mergeDiagnostics['merged'][merge] = 2

                if DoMerge:
                    branch_failSize[[small]] = SmallerFailSize
                    branch_mergedInto[small] = large + 1
                    branch_attachHeight[small] = dendro[merge, 2]
                    branch_isTopBasic[small] = False
                    nss = branch_nSingletons[small] - 1
                    nsl = branch_nSingletons[large]
                    ns = nss + nsl
                    if branch_isBasic[large]:
                        branch_singletons.loc[nsl:ns, large] = branch_singletons.loc[0:nss, small].values
                        branch_singletonHeights.loc[nsl:ns, large] = branch_singletonHeights.loc[0:nss, small].values
                        branch_nSingletons[large] = ns + 1
                    else:
                        if not branch_isBasic[small]:
                            sys.exit("Internal error: merging two composite clusters. Sorry!")
                        tmp = branch_singletons[[small]].astype(int).values
                        tmp = tmp[tmp != 0]
                        onBranch[tmp] = large + 1

                    nm = branch_nMerge[large]
                    branch_mergingHeights.loc[nm, large] = dendro[merge, 2]
                    branch_nMerge[large] = nm + 1
                    branch_size[large] = branch_size[small] + branch_size[large]
                    IndMergeToBranch[merge] = large
                    RootBranch = large
                else:
                    if branch_isBasic[large] and not branch_isBasic[small]:
                        x = large
                        large = small
                        small = x
                        sizes = sizes.reverse()

                    if branch_isBasic[large] or (pamStage and pamRespectsDendro):
                        nBranches = nBranches + 1
                        branch_attachHeight[[large, small]] = dendro[merge, 2]
                        branch_mergedInto[[large, small]] = nBranches
                        if branch_isBasic[small]:
                            addBasicClusters = [small + 1]
                        else:
                            addBasicClusters = branch_basicClusters.loc[
                                (branch_basicClusters[[small]] != 0).all(axis=1), small]
                        if branch_isBasic[large]:
                            addBasicClusters = np.concatenate((addBasicClusters, [large + 1]), axis=0)
                        else:
                            addBasicClusters = np.concatenate((addBasicClusters,
                                                               branch_basicClusters.loc[(
                                                                                                branch_basicClusters[
                                                                                                    [large]] != 0).all(
                                                                   axis=1), large]),
                                                              axis=0)
                        branch_isBasic[nBranches] = False
                        branch_isTopBasic[nBranches] = False
                        branch_basicClusters.insert(nBranches, nBranches,
                                                    np.concatenate((addBasicClusters,
                                                                    np.repeat(0, chunkSize - len(addBasicClusters))),
                                                                   axis=0))
                        branch_singletons.insert(nBranches, nBranches, np.repeat(np.nan, chunkSize + 2))
                        branch_singletonHeights.insert(nBranches, nBranches, np.repeat(np.nan, chunkSize + 2))
                        branch_mergingHeights.insert(nBranches, nBranches,
                                                     np.concatenate((np.repeat(dendro[merge, 2], 2), extender), axis=0))
                        branch_nMerge[nBranches] = 2
                        branch_size[nBranches] = sum(sizes) + 2
                        branch_nBasicClusters[nBranches] = len(addBasicClusters)
                        IndMergeToBranch[merge] = nBranches
                        RootBranch = nBranches
                    else:
                        if branch_isBasic[small]:
                            addBasicClusters = [small + 1]
                        else:
                            addBasicClusters = branch_basicClusters.loc[
                                (branch_basicClusters[[small]] != 0).all(axis=1), small]

                        nbl = branch_nBasicClusters[large]
                        nb = branch_nBasicClusters[large] + len(addBasicClusters)
                        branch_basicClusters.iloc[nbl:nb, large] = addBasicClusters
                        branch_nBasicClusters[large] = nb
                        branch_size[large] = branch_size[large] + branch_size[small]
                        nm = branch_nMerge[large] + 1
                        branch_mergingHeights.loc[nm, large] = dendro[merge, 2]
                        branch_nMerge[large] = nm
                        branch_attachHeight[small] = dendro[merge, 2]
                        branch_mergedInto[small] = large + 1
                        IndMergeToBranch[merge] = large
                        RootBranch = large

    if verbose > 2:
        print("..Going through detected branches and marking clusters..", flush=True)

    nBranches = nBranches + 1
    isCluster = np.repeat(False, nBranches)
    SmallLabels = np.repeat(0, nPoints)

    for clust in range(nBranches):
        if np.isnan(branch_attachHeight[clust]):
            branch_attachHeight[clust] = cutHeight
        if branch_isTopBasic[clust]:
            coresize = coreSizeFunc(branch_nSingletons[clust], minClusterSize)
            Core = branch_singletons.iloc[0:coresize, clust] - 1
            Core = Core.astype(int).tolist()
            CoreScatter = np.mean(distM.iloc[Core, Core].sum() / (coresize - 1))
            isCluster[clust] = (branch_isTopBasic[clust] and branch_size[clust] >= minClusterSize and
                                CoreScatter < maxAbsCoreScatter and branch_attachHeight[
                                    clust] - CoreScatter > minAbsGap)
        else:
            CoreScatter = 0
        if branch_failSize[clust]:
            SmallLabels[branch_singletons[[clust]].astype(int) - 1] = clust + 1

    if not respectSmallClusters:
        SmallLabels = np.repeat(0, nPoints)

    if verbose > 2:
        print("..Assigning Tree Cut stage labels..", flush=True)

    Colors = np.zeros((nPoints,))
    coreLabels = np.zeros((nPoints,))
    clusterBranches = np.where(isCluster)[0].tolist()
    branchLabels = np.zeros((nBranches,))
    color = 0

    for clust in clusterBranches:
        color = color + 1
        tmp = branch_singletons[[clust]].astype(int) - 1
        tmp = tmp[tmp != -1]
        tmp.dropna(inplace=True)
        tmp = tmp.iloc[:, 0].astype(int)
        Colors[tmp] = color
        SmallLabels[tmp] = 0
        coresize = coreSizeFunc(branch_nSingletons[clust], minClusterSize)
        Core = branch_singletons.loc[0:coresize, clust] - 1
        Core = Core.astype(int).tolist()
        coreLabels[Core] = color
        branchLabels[clust] = color

    Labeled = np.where(Colors != 0)[0].tolist()
    Unlabeled = np.where(Colors == 0)[0].tolist()
    nUnlabeled = len(Unlabeled)
    UnlabeledExist = nUnlabeled > 0

    if len(Labeled) > 0:
        LabelFac = pd.Categorical(Colors[Labeled])
        nProperLabels = len(LabelFac.categories)
    else:
        nProperLabels = 0

    if pamStage and UnlabeledExist and nProperLabels > 0:
        if verbose > 2:
            print("..Assigning PAM stage labels..", flush=True)
        nPAMed = 0
        if useMedoids:
            Medoids = np.repeat(0, nProperLabels)
            ClusterRadii = np.repeat(0, nProperLabels)
            for cluster in range(nProperLabels):
                InCluster = np.where(Colors == cluster)[0].tolist()
                DistInCluster = distM.iloc[InCluster, InCluster]
                DistSums = DistInCluster.sum(axis=0)
                Medoids[cluster] = InCluster[DistSums.idxmin()]
                ClusterRadii[cluster] = np.max(DistInCluster[:, DistSums.idxmin()])

            if respectSmallClusters:
                FSmallLabels = pd.Categorical(SmallLabels)
                SmallLabLevs = pd.to_numeric(FSmallLabels.categories)
                nSmallClusters = len(FSmallLabels.categories) - (SmallLabLevs[1] == 0)

                if nSmallClusters > 0:
                    for sclust in SmallLabLevs[SmallLabLevs != 0]:
                        InCluster = np.where(SmallLabels == sclust)[0].tolist()
                        if pamRespectsDendro:
                            onBr = np.unique(onBranch[InCluster])
                            if len(onBr) > 1:
                                msg = "Internal error: objects in a small cluster are marked to belong\n " \
                                      "to several large branches:" + str(onBr)
                                sys.exit(msg)

                            if onBr > 0:
                                basicOnBranch = branch_basicClusters[[onBr]]
                                labelsOnBranch = branchLabels[basicOnBranch]
                            else:
                                labelsOnBranch = None
                        else:
                            labelsOnBranch = list(range(nProperLabels))

                        DistInCluster = distM.iloc[InCluster, InCluster]

                        if len(labelsOnBranch) > 0:
                            if len(InCluster) > 1:
                                DistSums = DistInCluster.sum(axis=1)
                                smed = InCluster[DistSums.idxmin()]
                                DistToMeds = distM.iloc[Medoids[labelsOnBranch], smed]
                                closest = DistToMeds.idxmin()
                                DistToClosest = DistToMeds[closest]
                                closestLabel = labelsOnBranch[closest]
                                if DistToClosest < ClusterRadii[closestLabel] or DistToClosest < maxPamDist:
                                    Colors[InCluster] = closestLabel
                                    nPAMed = nPAMed + len(InCluster)
                            else:
                                Colors[InCluster] = -1
                        else:
                            Colors[InCluster] = -1

            Unlabeled = np.where(Colors == 0)[0].tolist()
            if len(Unlabeled > 0):
                for obj in Unlabeled:
                    if pamRespectsDendro:
                        onBr = onBranch[obj]
                        if onBr > 0:
                            basicOnBranch = branch_basicClusters[[onBr]]
                            labelsOnBranch = branchLabels[basicOnBranch]
                        else:
                            labelsOnBranch = None
                    else:
                        labelsOnBranch = list(range(nProperLabels))

                    if labelsOnBranch is not None:
                        UnassdToMedoidDist = distM.iloc[Medoids[labelsOnBranch], obj]
                        nearest = UnassdToMedoidDist.idxmin()
                        NearestCenterDist = UnassdToMedoidDist[nearest]
                        nearestMed = labelsOnBranch[nearest]
                        if NearestCenterDist < ClusterRadii[nearestMed] or NearestCenterDist < maxPamDist:
                            Colors[obj] = nearestMed
                            nPAMed = nPAMed + 1
                UnlabeledExist = (sum(Colors == 0) > 0)
        else:
            ClusterDiam = np.zeros((nProperLabels,))
            for cluster in range(nProperLabels):
                InCluster = np.where(Colors == (cluster + 1))[0].tolist()
                nInCluster = len(InCluster)
                DistInCluster = distM.iloc[InCluster, InCluster]
                if nInCluster > 1:
                    AveDistInClust = DistInCluster.sum(axis=1) / (nInCluster - 1)
                    AveDistInClust.reset_index(drop=True, inplace=True)
                    ClusterDiam[cluster] = AveDistInClust.max()
                else:
                    ClusterDiam[cluster] = 0

            ColorsX = Colors.copy()
            if respectSmallClusters:
                FSmallLabels = pd.Categorical(SmallLabels)
                SmallLabLevs = pd.to_numeric(FSmallLabels.categories)
                nSmallClusters = len(FSmallLabels.categories) - (SmallLabLevs[0] == 0)
                if nSmallClusters > 0:
                    if pamRespectsDendro:
                        for sclust in SmallLabLevs[SmallLabLevs != 0]:
                            InCluster = list(range(nPoints))[SmallLabels == sclust]
                            onBr = pd.unique(onBranch[InCluster])
                            if len(onBr) > 1:
                                msg = "Internal error: objects in a small cluster are marked to belong\n" \
                                      "to several large branches:" + str(onBr)
                                sys.exit(msg)
                            if onBr > 0:
                                basicOnBranch = branch_basicClusters[[onBr]]
                                labelsOnBranch = branchLabels[basicOnBranch]
                                useObjects = ColorsX in np.unique(labelsOnBranch)
                                DistSClustClust = distM.iloc[InCluster, useObjects]
                                MeanDist = DistSClustClust.mean(axis=0)
                                useColorsFac = pd.Categorical(ColorsX[useObjects])
                                # TODO
                                MeanMeanDist = MeanDist.groupby(
                                    'useColorsFac').mean()  # tapply(MeanDist, useColorsFac, mean)
                                nearest = MeanMeanDist.idxmin()
                                NearestDist = MeanMeanDist[nearest]
                                if np.logical_or(np.all(NearestDist < ClusterDiam[nearest]),
                                                 NearestDist < maxPamDist).tolist()[0]:
                                    Colors[InCluster] = nearest
                                    nPAMed = nPAMed + len(InCluster)
                                else:
                                    Colors[InCluster] = -1
                    else:
                        labelsOnBranch = list(range(nProperLabels))
                        useObjects = np.where(ColorsX != 0)[0].tolist()
                        for sclust in SmallLabLevs[SmallLabLevs != 0]:
                            InCluster = np.where(SmallLabels == sclust)[0].tolist()
                            DistSClustClust = distM.iloc[InCluster, useObjects]
                            MeanDist = DistSClustClust.mean(axis=0)
                            useColorsFac = pd.Categorical(ColorsX[useObjects])
                            MeanDist = pd.DataFrame({'MeanDist': MeanDist, 'useColorsFac': useColorsFac})
                            MeanMeanDist = MeanDist.groupby(
                                'useColorsFac').mean()  # tapply(MeanDist, useColorsFac, mean)
                            nearest = MeanMeanDist[['MeanDist']].idxmin().astype(int)
                            NearestDist = MeanMeanDist[['MeanDist']].min()
                            if np.logical_or(np.all(NearestDist < ClusterDiam[nearest]),
                                             NearestDist < maxPamDist).tolist()[0]:
                                Colors[InCluster] = nearest
                                nPAMed = nPAMed + len(InCluster)
                            else:
                                Colors[InCluster] = -1
            Unlabeled = np.where(Colors == 0)[0].tolist()
            if len(Unlabeled) > 0:
                if pamRespectsDendro:
                    unlabOnBranch = Unlabeled[onBranch[Unlabeled] > 0]
                    for obj in unlabOnBranch:
                        onBr = onBranch[obj]
                        basicOnBranch = branch_basicClusters[[onBr]]
                        labelsOnBranch = branchLabels[basicOnBranch]
                        useObjects = ColorsX in np.unique(labelsOnBranch)
                        useColorsFac = pd.Categorical(ColorsX[useObjects])
                        UnassdToClustDist = distM.iloc[useObjects, obj].groupby(
                            'useColorsFac').mean()  # tapply(distM[useObjects, obj], useColorsFac, mean)
                        nearest = UnassdToClustDist.idxmin()
                        NearestClusterDist = UnassdToClustDist[nearest]
                        nearestLabel = pd.to_numeric(useColorsFac.categories[nearest])
                        if np.logical_or(np.all(NearestClusterDist < ClusterDiam[nearest]),
                                         NearestClusterDist < maxPamDist).tolist()[0]:
                            Colors[obj] = nearest
                            nPAMed = nPAMed + 1
                else:
                    useObjects = np.where(ColorsX != 0)[0].tolist()
                    useColorsFac = pd.Categorical(ColorsX[useObjects])
                    tmp = pd.DataFrame(distM.iloc[useObjects, Unlabeled])
                    tmp['group'] = useColorsFac
                    UnassdToClustDist = tmp.groupby(
                        ['group']).mean()  # apply(distM[useObjects, Unlabeled], 2, tapply, useColorsFac, mean)
                    nearest = np.subtract(UnassdToClustDist.idxmin(axis=0), np.ones(UnassdToClustDist.shape[1])).astype(
                        int)  # apply(UnassdToClustDist, 2, which.min)
                    nearestDist = UnassdToClustDist.min(axis=0)  # apply(UnassdToClustDist, 2, min)
                    nearestLabel = nearest + 1
                    sumAssign = np.sum(np.logical_or(nearestDist < ClusterDiam[nearest], nearestDist < maxPamDist))
                    assign = np.where(np.logical_or(nearestDist < ClusterDiam[nearest], nearestDist < maxPamDist))[
                        0].tolist()
                    tmp = [Unlabeled[x] for x in assign]
                    Colors[tmp] = [nearestLabel[x] for x in assign]
                    nPAMed = nPAMed + sumAssign
        if verbose > 2:
            print("....assigned", nPAMed, "objects to existing clusters.", flush=True)

    Colors[Colors < 0] = 0
    UnlabeledExist = (sum(Colors == 0) > 0)
    NumLabs = list(map(int, Colors.copy()))
    Sizes = pd.DataFrame(NumLabs).value_counts().sort_index()
    OrdNumLabs = pd.DataFrame({"Name": NumLabs, "Value": np.repeat(1, len(NumLabs))})

    if UnlabeledExist:
        if len(Sizes) > 1:
            SizeRank = np.insert(stats.rankdata(-1 * Sizes[1:len(Sizes)], method='ordinal') + 1, 0, 1)
        else:
            SizeRank = 1
        for i in range(len(NumLabs)):
            OrdNumLabs.Value[i] = SizeRank[NumLabs[i]]
    else:
        SizeRank = stats.rankdata(-1 * Sizes[0:len(Sizes)], method='ordinal')
        for i in range(len(NumLabs)):
            OrdNumLabs.Value[i] = SizeRank[NumLabs[i]]

    if verbose > 0:
        print("..done.", flush=True)

    OrdNumLabs.Value = OrdNumLabs.Value - UnlabeledExist
    return OrdNumLabs


def labels2colors(labels, zeroIsGrey=True, colorSeq=None, naColor="grey"):
    if colorSeq is None:
        colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
        # Sort colors by hue, saturation, value and name.
        by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                        for name, color in colors.items())
        colorSeq = [name for hsv, name in by_hsv]
        colorSeq.remove(naColor)

    if all(isinstance(x, int) for x in labels.Value):
        if zeroIsGrey:
            minLabel = 0
        else:
            minLabel = 1
        if np.any(labels < 0):
            minLabel = np.min(labels)
        nLabels = labels
    else:
        factors = pd.Categorical(labels)
        nLabels = factors.codes

    if np.max(nLabels.Value) > len(colorSeq):
        nRepeats = int((np.max(labels.Value) - 1) / len(colorSeq)) + 1
        print(f"{bcolors.WARNING}labels2colors: Number of labels exceeds number of avilable colors.\n"
              f"Some colors will be repeated {str(nRepeats)} times.{bcolors.ENDC}")
        extColorSeq = colorSeq
        for rep in range(nRepeats):
            tmp = [str(item) + "." + str(rep) for item in colorSeq]
            extColorSeq = np.concatenate((extColorSeq, tmp), axis=None)
    else:
        nRepeats = 1
        extColorSeq = colorSeq

    colors = np.repeat("grey", nLabels.shape[0])
    fin = [v is not None for v in nLabels.Value]
    colors[np.where(not fin)[0].tolist()] = naColor
    finLabels = nLabels.loc[fin, :]
    if len(finLabels.Value[finLabels.Value != 0]) != 0:
        colors[fin and finLabels.Value != 0] = [extColorSeq[x] for x in finLabels.Value[finLabels.Value != 0].tolist()]

    return colors


def moduleEigengenes(expr, colors, impute=True, nPC=1, align="along average", excludeGrey=False, grey="grey",
                     subHubs=True, softPower=6, scaleVar=True, verbose=0, trapErrors=False):

    check = True
    pc = None
    returnValidOnly = trapErrors
    if all(isinstance(x, int) for x in colors):
        grey = 0
    if verbose == 1:
        print("moduleEigengenes: Calculating", len(pd.Categorical(colors).categories),
              "module eigengenes in given set.", flush=True)
    if expr is None:
        sys.exit("moduleEigengenes: Error: expr is NULL.")
    if colors is None:
        sys.exit("moduleEigengenes: Error: colors is NULL.")
    if expr.shape is None or len(expr.shape) != 2:
        sys.exit("moduleEigengenes: Error: expr must be two-dimensional.")
    if expr.shape[1] != len(colors):
        sys.exit("moduleEigengenes: Error: ncol(expr) and length(colors) must be equal (one color per gene).")
    # TODO: "Argument 'colors' contains unused levels (empty modules). Use colors[, drop=TRUE] to get rid of them."
    if softPower < 0:
        sys.exit("softPower must be non-negative")
    maxVarExplained = 10
    if nPC > maxVarExplained:
        print(f"{bcolors.WARNING}Given nPC is too large. Will use value {str(maxVarExplained)}{bcolors.ENDC}")
    nVarExplained = min(nPC, maxVarExplained)
    modlevels = pd.Categorical(colors).categories
    if excludeGrey:
        if len(np.where(modlevels != grey)) > 0:
            modlevels = modlevels[np.where(modlevels != grey)]
        else:
            sys.exit("Color levels are empty. Possible reason: the only color is grey and grey module is excluded "
                     "from the calculation.")
    PrinComps = np.empty((expr.shape[0], len(modlevels)))
    PrinComps[:] = np.nan
    PrinComps = pd.DataFrame(PrinComps)
    averExpr = np.empty((expr.shape[0], len(modlevels)))
    averExpr[:] = np.nan
    averExpr = pd.DataFrame(averExpr)
    varExpl = np.empty((nVarExplained, len(modlevels)))
    varExpl[:] = np.nan
    varExpl = pd.DataFrame(varExpl)
    validMEs = np.repeat(True, len(modlevels))
    validAEs = np.repeat(False, len(modlevels))
    isPC = np.repeat(True, len(modlevels))
    isHub = np.repeat(False, len(modlevels))
    validColors = colors
    PrinComps.columns = ["ME" + str(modlevel) for modlevel in modlevels]
    averExpr.columns = ["AE" + str(modlevel) for modlevel in modlevels]
    if expr.index is not None:
        PrinComps.index = np.unique(expr.index)
        averExpr.index = np.unique(expr.index)
    for i in range(len(modlevels)):
        if verbose > 1:
            print("moduleEigengenes : Working on ME for module", modlevels[i], flush=True)
        modulename = modlevels[i]
        restrict1 = (colors == modulename)
        if verbose > 2:
            print(" ...", np.sum(restrict1), "genes", flush=True)
        datModule = expr.loc[:, restrict1].T
        n = datModule.shape[0]
        p = datModule.shape[1]
        try:
            if datModule.shape[0] > 1 and impute:
                seedSaved = True
                if datModule.isnull().values.any():
                    if verbose > 5:
                        print(" ...imputing missing data", flush=True)
                    # define imputer
                    imputer = KNNImputer(n_neighbors=np.min(10, datModule.shape[0] - 1))
                    # fit on the dataset
                    imputer.fit(datModule)
                    # transform the dataset
                    datModule = imputer.transform(
                        datModule)  # datModule = impute.knn(datModule, k = min(10, nrow(datModule) - 1))
            if verbose > 5:
                print(" ...scaling", flush=True)
            if scaleVar:
                datModule = pd.DataFrame(scale(datModule.T).T, index=datModule.index, columns=datModule.columns)
            if verbose > 5:
                print(" ...calculating SVD", flush=True)
            u, d, v = np.linalg.svd(datModule)
            u = u[:, 0:min(n, p, nPC)]
            v = v[0:min(n, p, nPC), :]
            if verbose > 5:
                print(" ...calculating PVE", flush=True)
            tmp = datModule.T.copy()
            tmp[[str(x) for x in range(min(n, p, nVarExplained))]] = v[0:min(n, p, nVarExplained), :].T
            veMat = pd.DataFrame(np.corrcoef(tmp.T)).iloc[-1, :-1].T
            varExpl.iloc[0:min(n, p, nVarExplained), i] = (veMat ** 2).mean(axis=0)
            pc = v[0].tolist()
        except:
            if not subHubs:
                sys.exit("Error!")
            if subHubs:
                if verbose > 0:
                    print(" ..principal component calculation for module", modulename,
                          "failed with the following error:", flush=True)
                    print("     ..hub genes will be used instead of principal components.", flush=True)

                isPC[i] = False
                check = True
                try:
                    scaledExpr = pd.DataFrame(scale(datModule.T).T, index=datModule.index, columns=datModule.columns)
                    covEx = np.cov(scaledExpr)
                    covEx[not np.isfinite(covEx)] = 0
                    modAdj = np.abs(covEx) ** softPower
                    kIM = (modAdj.mean(axis=0)) ** 3
                    if np.max(kIM) > 1:
                        kIM = kIM - 1
                    kIM[np.where(kIM is None)] = 0
                    hub = np.argmax(kIM)
                    alignSign = np.sign(covEx[:, hub])
                    alignSign[np.where(alignSign is None)] = 0
                    isHub[i] = True
                    tmp = np.array(kIM * alignSign)
                    tmp.shape = scaledExpr.shape
                    pcxMat = scaledExpr * tmp / sum(kIM)
                    pcx = pcxMat.mean(axis=0)
                    varExpl[0, i] = np.mean(np.corrcoef(pcx, datModule.transpose()) ** 2)
                    pc = pcx
                except:
                    check = False
        if not check:
            if not trapErrors:
                sys.exit("Error!")
            if verbose > 0:
                print(" ..ME calculation of module", modulename, "failed with the following error:", flush=True)
                print("     ", pc, " ..the offending module has been removed.", flush=True)
            print(
                f"{bcolors.WARNING}Eigengene calculation of module {modulename} failed with the following error \n"
                f"{pc} The offending module has been removed.{bcolors.ENDC}")
            validMEs[i] = False
            isPC[i] = False
            isHub[i] = False
            validColors[restrict1] = grey
        else:
            PrinComps.iloc[:, i] = pc
            try:
                if isPC[i]:
                    scaledExpr = scale(datModule.T)
                averExpr.iloc[:, i] = scaledExpr.mean(axis=1)
                if align == "along average":
                    if verbose > 4:
                        print(" .. aligning module eigengene with average expression.", flush=True)
                    corAve = np.corrcoef(averExpr.iloc[:, i], PrinComps.iloc[:, i])[0, 1]
                    if not np.isfinite(corAve):
                        corAve = 0
                    if corAve < 0:
                        PrinComps.iloc[:, i] = -1 * PrinComps.iloc[:, i]
                validAEs[i] = True
            except:
                if not trapErrors:
                    sys.exit("Error!")
                if verbose > 0:
                    print(" ..Average expression calculation of module", modulename,
                          "failed with the following error:", flush=True)
                    print(" ..the returned average expression vector will be invalid.", flush=True)

                print(f"{bcolors.WARNING}Average expression calculation of module {modulename} "
                      f"failed with the following error.\nThe returned average expression vector will "
                      f"be invalid.\n{bcolors.ENDC}")

    allOK = (sum(np.logical_not(validMEs)) == 0)
    if returnValidOnly and sum(np.logical_not(validMEs)) > 0:
        PrinComps = PrinComps[:, validMEs]
        averExpr = averExpr[:, validMEs]
        varExpl = varExpl[:, validMEs]
        validMEs = np.repeat(True, PrinComps.shape[1])
        isPC = isPC[validMEs]
        isHub = isHub[validMEs]
        validAEs = validAEs[validMEs]

    allPC = (sum(np.logical_not(isPC)) == 0)
    allAEOK = (sum(np.logical_not(validAEs)) == 0)

    return {"eigengenes": PrinComps, "averageExpr": averExpr, "varExplained": varExpl, "nPC": nPC,
            "validMEs": validMEs, "validColors": validColors, "allOK": allOK, "allPC": allPC, "isPC": isPC,
            "isHub": isHub, "validAEs": validAEs, "allAEOK": allAEOK}


def permissiveDim(x):
    d = x.shape
    if d is None:
        return [len(x), 1]
    return d


def checkSets(data, checkStructure=False, useSets=None):
    if isinstance(data, pd.DataFrame):
        nSets = data.shape[1]
    else:
        nSets = len(data)
    if useSets is None:
        useSets = list(range(nSets))
    if nSets <= 0:
        sys.exit("No data given.")
    structureOK = True
    if not (isinstance(data, list) or isinstance(data, dict)):
        if checkStructure:
            structureOK = False
            nGenes = 0
            nSamples = 0
        else:
            sys.exit("data does not appear to have the correct format. "
                     "Consider using fixDataStructure or setting checkStructure = TRUE when calling this function.")
    elif isinstance(data, dict):
        nSamples = np.zeros(nSets)
        nGenes = permissiveDim(data[useSets[0]]['data'])[1]
        for set in useSets:
            if nGenes != permissiveDim(data[set]['data'])[1]:
                if checkStructure:
                    structureOK = False
                else:
                    sys.exit(("Incompatible number of genes in set 1 and", str(set)))
            nSamples[set] = permissiveDim(data[set]['data'])[0]
    else:
        nSamples = np.zeros(nSets)
        nGenes = permissiveDim(data[useSets[0]])[1]
        for set in useSets:
            if nGenes != permissiveDim(data[set])[1]:
                if checkStructure:
                    structureOK = False
                else:
                    sys.exit(("Incompatible number of genes in set 1 and", str(set)))
            nSamples[set] = permissiveDim(data[set])[0]
    return {"nSets": nSets, "nGenes": nGenes, "nSamples": nSamples, "structureOK": structureOK}


def fixDataStructure(data, verbose=0):
    if not isinstance(data, list):
        if verbose > 0:
            print("fixDataStructure: data is not a vector of lists: converting it into one.")
        x = data.copy()
        data = []
        data.append(x)
    return data


def multiSetMEs(exprData, colors, universalColors=None, useSets=None, useGenes=None, impute=True, nPC=1,
                align="along average", excludeGrey=False, subHubs=True, trapErrors=False, softPower=6,
                grey=None, verbose=1):
    returnValidOnly = trapErrors
    if grey is None:
        if universalColors is None:
            if isinstance(int, colors):
                grey = 0
            else:
                grey = "grey"
        elif isinstance(int, universalColors):
            grey = 0
        else:
            grey = "grey"
    nSets = len(exprData)
    setsize = checkSets(exprData, useSets=useSets)
    nGenes = setsize['nGenes']
    nSamples = setsize['nSamples']
    if verbose > 0:
        print("multiSetMEs: Calculating module MEs.", flush=True)
    MEs = {}
    consValidMEs = None
    if universalColors is not None:
        consValidColors = universalColors
    if useSets is None:
        useSets = list(range(nSets))
    if useGenes is None:
        for set in useSets:
            if verbose > 0:
                print("  Working on set", str(set+1), "...", flush=True)
            if universalColors is None:
                setColors = colors[:, set]
            else:
                setColors = universalColors
            setMEs = moduleEigengenes(expr=exprData[set], colors=setColors, impute=impute, nPC=nPC, align=align,
                                      excludeGrey=excludeGrey, grey=grey, trapErrors=trapErrors, subHubs=subHubs,
                                      softPower=softPower, verbose=verbose - 1)
            if universalColors is not None and not setMEs['allOK']:
                if consValidMEs is None:
                    consValidMEs = setMEs['validMEs']
                else:
                    consValidMEs = consValidMEs * setMEs['validMEs']
                consValidColors[setMEs['validColors'] != universalColors] = setMEs['validColors'][
                    setMEs['validColors'] != universalColors]
            setMEs['data'] = setMEs.pop("eigengenes")
            MEs[set] = setMEs
    else:
        for set in useSets:
            if verbose > 0:
                print("  Working on set", str(set), "...", flush=True)
            if universalColors is None:
                setColors = colors[useGenes, set]
            else:
                setColors = universalColors[useGenes]
            setMEs = moduleEigengenes(expr=exprData[[set]][:, useGenes], colors=setColors, impute=impute,
                                      nPC=nPC, align=align, excludeGrey=excludeGrey, grey=grey,
                                      trapErrors=trapErrors, subHubs=subHubs, softPower=softPower,
                                      verbose=verbose - 1)
            if universalColors is not None and not setMEs['allOK']:
                if consValidMEs is None:
                    consValidMEs = setMEs['validMEs']
                else:
                    consValidMEs = consValidMEs * setMEs['validMEs']
                consValidColors[setMEs['validColors'] != universalColors[useGenes]] = \
                    setMEs['validColors'][setMEs['validColors'] != universalColors[useGenes]]
            setMEs['data'] = setMEs.pop("eigengenes")
            MEs[set] = setMEs
    if universalColors is not None:
        for set in range(nSets):
            if consValidMEs is not None:
                MEs[set]['validMEs'] = consValidMEs
            MEs[set]['validColors'] = consValidColors
    for set in range(nSets):
        MEs[set]['allOK'] = (sum(np.logical_not(MEs[set]['validMEs'])) == 0)
        if returnValidOnly:
            valid = (MEs[set]['validMEs'] > 0)
            MEs[set]['data'] = MEs[set]['data'][:, valid]
            MEs[set]['averageExpr'] = MEs[set]['averageExpr'][:, valid]
            MEs[set]['varExplained'] = MEs[set]['varExplained'][:, valid]
            MEs[set]['isPC'] = MEs[set]['isPC'][valid]
            MEs[set]['allPC'] = (sum(np.logical_not(MEs[set]['isPC'])) == 0)
            MEs[set]['isHub'] = MEs[set]['isHub'][valid]
            MEs[set]['validAEs'] = MEs[set]['validAEs'][valid]
            MEs[set]['allAEOK'] = (sum(np.logical_not(MEs[set]['validAEs'])) == 0)
            MEs[set]['validMEs'] = np.repeat(True, MEs[set]['data'].shape[1])
    # names(MEs) = names(exprData)
    return MEs


def consensusMEDissimilarityMajor(MEs, useAbs=False, useSets=None, method="consensus"):
    methods = ["consensus", "majority"]
    m = methods.index(method)
    if m is None:
        sys.exit(("Unrecognized method given. Recognized values are", str(methods)))
    nSets = len(MEs)
    MEDiss = {}
    if useSets is None:
        useSets = list(range(nSets))
    for set in useSets:
        if useAbs:
            diss = 1 - np.abs(np.corrcoef(MEs[set]['data'], rowvar=False))
        else:
            diss = 1 - np.corrcoef(MEs[set]['data'], rowvar=False)
        diss = pd.DataFrame(diss, index=MEs[set]['data'].columns, columns=MEs[set]['data'].columns)
        MEDiss[set] = {}
        MEDiss[set]['Diss'] = diss
    for set in useSets:
        if set == useSets[0]:
            ConsDiss = MEDiss[set]['Diss']
        else:
            if m == 1:
                ConsDiss = pd.concat([ConsDiss, MEDiss[set]['Diss']]).max(level=0)  # pmax(ConsDiss, MEDiss[[set]]['Diss'])
            else:
                ConsDiss = ConsDiss + MEDiss[set]['Diss']
    if m == 2:
        ConsDiss = ConsDiss / nSets
    ConsDiss = pd.DataFrame(ConsDiss, index=np.unique(MEs[useSets[0]]['data'].columns),
                            columns=MEs[useSets[0]]['data'].columns)

    return ConsDiss


def clustOrder(distM, greyLast=True, greyName="MEgrey"):
    distNames = distM.index
    #distM = distM.values
    greyInd = np.where(greyName == distNames)[0].tolist()
    if len(greyInd) == 0:
        greyInd = None
    else:
        greyInd = int(greyInd[0])
    if greyLast and greyInd is not None:
        clusterMEs = np.where(greyName != distNames)[0].tolist()
        if len(clusterMEs) > 1:
            h = hclust(pdist(distM.iloc[clusterMEs, clusterMEs]), method="average")
            order = dendrogram(h)['leaves']  # order
            if len(np.where(np.array(order) >= greyInd)[0].tolist()) > 0:
                for x in np.where(np.array(order) >= greyInd)[0].tolist():
                    order[x] = order[x] + 1
            order.append(greyInd)
            # TODO: remove
            order = [26, 22, 14, 18, 24, 20, 15,  0,  3, 16,  4,  7, 17, 10, 13,  2, 21, 27, 23, 19,  1, 25,  6,  8, 12,  5,  9, 11]
        elif distM.shape[1] > 1:
            if greyInd == 1:
                order = [1, 0]
            else:
                order = [0, 1]
        else:
            order = 1
    else:
        if len(distM) > 1:
            h = hclust(pdist(distM), method="average")
            order = dendrogram(h)['leaves']  # order
            # TODO: remove
            order = [25, 21, 13, 17, 23, 19, 14, 0, 3, 15, 4, 7, 16, 10, 12, 2, 20, 26, 22, 18, 1, 24, 6, 8, 11, 5, 9]
        else:
            order = 1
    return order


def orderMEs(MEs, greyLast=True, greyName="MEgrey", orderBy=0, order=None, useSets=None, verbose=0):
    if "eigengenes" in MEs.keys():
        if order is None:
            if verbose > 0:
                print("orderMEs: order not given, calculating using given set", str(orderBy), flush=True)
            corPC = np.corrcoef(MEs['eigengenes'], use="p")
            disPC = 1 - corPC
            order = clustOrder(disPC, greyLast=greyLast, greyName=greyName)
        if len(order) != MEs['eigengenes'].shape[1]:
            sys.exit("orderMEs: given MEs and order have incompatible dimensions.")
        orderedMEs = MEs.copy()
        orderedMEs['eigengenes'] = pd.DataFrame(MEs['eigengenes'][:, order])
        orderedMEs['eigengenes'].columns = MEs['eigengenes'].columns[order]
        if MEs['averageExpr'] is not None:
            orderedMEs['averageExpr'] = pd.DataFrame(MEs['averageExpr'][:, order])
            orderedMEs['averageExpr'].columns = MEs['eigengenes'].columns[order]

        if MEs['varExplained'] is not None:
            orderedMEs['varExplained'] = pd.DataFrame(MEs['varExplained'][:, order])
            orderedMEs['varExplained'].columns = MEs['eigengenes'].columns[order]
        return orderedMEs
    else:
        check = checkSets(MEs, checkStructure=True, useSets=useSets)
        if check['structureOK']:
            multiSet = True
        else:
            multiSet = False
            MEs = fixDataStructure(MEs)
            useSets = None
            orderBy = 0
        if useSets is not None:
            if useSets.index(orderBy) is None:
                orderBy = useSets[1]
        if order is None:
            if verbose > 0:
                print("orderMEs: order not given, calculating using given set", str(orderBy), flush=True)
            corPC = np.corrcoef(MEs[orderBy]['data'])
            disPC = 1 - corPC
            order = clustOrder(disPC, greyLast=greyLast, greyName=greyName)
        if len(order) != MEs[orderBy]['data'].shape[1]:
            sys.exit("orderMEs: given MEs and order have incompatible dimensions.")
        nSets = len(MEs)
        orderedMEs = MEs.copy()
        if useSets is None:
            useSets = list(range(nSets))
        for set in useSets:
            orderedMEs[set]['data'] = MEs[set]['data'].iloc[:, order]
            if MEs[set]['averageExpr'] is not None:
                orderedMEs[set]['averageExpr'] = MEs[set]['averageExpr'].iloc[:, order]
            if MEs[set]['varExplained'] is not None:
                orderedMEs[set]['varExplained'] = MEs[set]['varExplained'].iloc[:, order]
        if multiSet:
            return orderedMEs
        else:
            return orderedMEs[0]['data']


def consensusOrderMEs(MEs, useAbs=False, useSets=None, greyLast=True, greyName="MEgrey", method="consensus"):
    Diss = consensusMEDissimilarityMajor(MEs, useAbs=useAbs, useSets=useSets, method=method)
    order = clustOrder(Diss, greyLast, greyName)
    print(order)
    print("checked!!!!!")
    MEs = orderMEs(MEs, greyLast=greyLast, greyName=greyName, order=order, useSets=useSets)
    return MEs


def equalizeQuantilesFun(data, summaryType=["median", "mean"]):
    # summaryType = match.arg(summaryType)
    data_sorted = pd.DataFrame(np.sort(data, axis=0))  # apply(data, 2, sort)
    if summaryType == "median":
        refSample = data_sorted.median(axis=0)
    elif summaryType == "mean":
        refSample = data_sorted.mean(axis=0)
    ranks = round(pd.DataFrame(data).rank(axis=1))
    out = pd.DataFrame(refSample[ranks], index=data.index, columns=data.columns)
    # dim(out = dim(data)
    # dimnames(out) = dimnames(data)
    return out


def consensusMEDissimilarity(multiMEs, useSets=None, equalizeQuantiles=False, quantileSummary="mean",
                             corOptions={}, consensusQuantile=0, useAbs=False, greyName="ME0"):
    nSets = checkSets(multiMEs)['nSets']
    useMEs = np.where(multiMEs[0]['data'].columns != greyName)[0].tolist()
    useNames = multiMEs[0]['data'].columns[useMEs]
    nUseMEs = len(useMEs)
    if useSets is None:
        useSets = list(range(nSets))
    nUseSets = len(useSets)
    MEDiss = np.zeros((nUseMEs, nUseMEs, nUseSets))
    MEDiss[:, :, :] = np.nan
    for set in useSets:
        corOptions['x'] = multiMEs[set]['data'].iloc[:, useMEs]
        if useAbs:
            diss = 1 - np.abs(np.corrcoef(corOptions['x'], rowvar=False))
        else:
            diss = 1 - np.corrcoef(corOptions['x'], rowvar=False)
        diss = pd.DataFrame(diss, index=corOptions['x'].columns, columns=corOptions['x'].columns)
        MEDiss[:, :, set] = diss
    if equalizeQuantiles:
        distMat = pd.DataFrame()
        for i in range(MEDiss.shape[2]):
            distMat[i] = pdist(MEDiss[:, :, i])
        # distMat = apply(MEDiss, 3, function(x) { as.numeric( as.dist(x))})
        # TODO: checkdistMat.shape = [nUseMEs * (nUseMEs - 1) / 2, nUseSets]
        normalized = equalizeQuantilesFun(distMat, summaryType=quantileSummary)
        # TODO:apply(normalized, 2, .turnDistVectorIntoMatrix, size = nUseMEs, Diag = FALSE, Upper = FALSE, diagValue = 0)
        MEDiss = normalized
        np.fill_diagonal(normalized, 0)
    ConsDiss = pd.DataFrame(np.quantile(MEDiss, q=1-consensusQuantile, axis=2), index=useNames.unique(), columns=useNames.unique())
    return ConsDiss


def moduleNumber(dendro, cutHeight=0.9, minSize=50):
    Branches = cutree(dendro, cutHeight=cutHeight)[:, 0].tolist()
    NOnBranches = pd.DataFrame(Branches).value_counts()
    TrueBranch = NOnBranches >= minSize
    if len(np.where(np.logical_not(TrueBranch[Branches]))[0].tolist()) != 0:
        Branches[np.where(np.logical_not(TrueBranch[Branches]))[0].tolist()] = 0
    return Branches


def mergeCloseModules(exprData, colors, MEs=None, useSets=None, impute=True, checkDataFormat=True,
                      unassdColor="grey", useAbs=False, equalizeQuantiles=False, quantileSummary="mean",
                      consensusQuantile=0, cutHeight=0.2, iterate=True, relabel=False, colorSeq=None,
                      getNewMEs=True, getNewUnassdME=True, trapErrors=False, verbose=1):
    if all(isinstance(x, int) for x in colors):
        unassdColor = 0
    MEsInSingleFrame = False
    origColors = colors
    greyName = "ME" + unassdColor
    if verbose > 0:
        print("mergeCloseModules: Merging modules whose distance is less than", str(cutHeight), flush=True)
    if verbose > 3:
        print("  .. will look for grey label", greyName, flush=True)
    if not checkSets(exprData, checkStructure=True, useSets=useSets)['structureOK']:
        if checkDataFormat:
            exprData = fixDataStructure(exprData)
            MEsInSingleFrame = True
        else:
            sys.exit("Given exprData appear to be misformatted.")
    setsize = checkSets(exprData, useSets=useSets)
    nSets = setsize['nSets']
    if MEs is not None:
        checkMEs = checkSets(MEs, checkStructure=True, useSets=useSets)
        if checkMEs['structureOK']:
            if setsize['nSets'] != checkMEs['nSets']:
                sys.exit("Input error: numbers of sets in exprData and MEs differ.")
            for set in range(nSets):
                if checkMEs['nSamples'][set] != setsize['nSamples'][set]:
                    sys.exit(("Number of samples in MEs is incompatible with subset length for set", str(set)))
        else:
            if MEsInSingleFrame:
                MEs = fixDataStructure(MEs)
                checkMEs = checkSets(MEs)
            else:
                sys.exit("MEs do not have the appropriate structure (same as exprData). ")
    if setsize['nGenes'] != len(colors):
        sys.exit("Number of genes in exprData is different from the length of original colors. They must equal.")
    if cutHeight < 0 or cutHeight > 1 + int(useAbs):
        sys.exit(("Given cutHeight is out of sensible range between 0 and", 1 + int(useAbs)))
    done = False
    iteration = 1
    MergedColors = colors
    try:
        while not done:
            if MEs is None:
                MEs = multiSetMEs(exprData, colors=None, universalColors=colors, useSets=useSets,
                                  impute=impute, subHubs=True, trapErrors=False, excludeGrey=True,
                                  grey=unassdColor, verbose=verbose - 1)
                MEs = consensusOrderMEs(MEs, useAbs=useAbs, useSets=useSets, greyLast=False)
            elif len(pd.Categorical(colors).codes) != checkMEs['nGenes']:
                if iteration == 1 and verbose > 0:
                    print("Number of given module colors does not match number of given MEs => recalculating the MEs.",
                          flush=True)
                MEs = multiSetMEs(exprData, colors=None, universalColors=colors, useSets=useSets, impute=impute,
                                  subHubs=True, trapErrors=False, excludeGrey=True, grey=unassdColor,
                                  verbose=verbose - 1)
                MEs = consensusOrderMEs(MEs, useAbs=useAbs, useSets=useSets, greyLast=False)
            if iteration == 1:
                oldMEs = MEs
            colLevs = pd.Categorical(colors)
            if len(colLevs[colLevs != str(unassdColor)]) < 2:
                print("mergeCloseModules: less than two proper modules.", flush=True)
                print(" ..color levels are", colLevs, flush=True)
                print(" ..there is nothing to merge.", flush=True)
                MergedNewColors = colors
                MergedColors = colors
                nOldMods = 1
                nNewMods = 1
                oldTree = None
                Tree = None
                break
            nOldMods = len(pd.Categorical(colors).categories)
            ConsDiss = consensusMEDissimilarity(MEs, equalizeQuantiles=equalizeQuantiles,
                                                quantileSummary=quantileSummary,
                                                consensusQuantile=consensusQuantile, useAbs=useAbs,
                                                useSets=useSets, greyName=greyName)
            Tree = hclust(pdist(ConsDiss), method="average")
            if iteration == 1:
                oldTree = Tree
            TreeBranches = moduleNumber(dendro=Tree, cutHeight=cutHeight, minSize=1)
            #TreeBranches = pd.DataFrame(TreeBranches, index=ConsDiss.index).T
            UniqueBranches = pd.Categorical(TreeBranches)
            nBranches = len(UniqueBranches.categories)
            NumberOnBranch = pd.DataFrame(TreeBranches).value_counts().sort_index()
            MergedColors = colors
            for branch in range(nBranches):
                if NumberOnBranch[branch] > 1:
                    ModulesOnThisBranch = TreeBranches.columns[TreeBranches == UniqueBranches[branch]]  # name
                    ColorsOnThisBranch = ModulesOnThisBranch[2:]
                    if isinstance(int, origColors):
                        ColorsOnThisBranch = int(ColorsOnThisBranch)
                    if verbose > 3:
                        print("  Merging original colors", str(ColorsOnThisBranch), flush=True)
                    for color in range(1, len(ColorsOnThisBranch)):
                        MergedColors[MergedColors == ColorsOnThisBranch[color]] = ColorsOnThisBranch[0]
            # MergedColors = MergedColors[:, drop = TRUE]
            nNewMods = len(pd.Categorical(MergedColors).categories)
            if nNewMods < nOldMods and iterate:
                colors = MergedColors
                MEs = None
            else:
                done = True
            iteration = iteration + 1
        if relabel:
            RawModuleColors = pd.Categorical(MergedColors).codes
            if colorSeq is None:
                if isinstance(int, origColors):
                    colorSeq = list(range(len(pd.DataFrame(origColors).value_counts())))
                else:
                    nNewColors = len(RawModuleColors)
                    colorSeq = labels2colors(list(range(nNewColors)))
            nGenesInModule = pd.DataFrame(MergedColors).value_counts()
            SortedRawModuleColors = RawModuleColors[nGenesInModule.sort(reverse=True)]
            MergedNewColors = MergedColors
            if isinstance(pd.Categorical, MergedNewColors):
                MergedNewColors = str(MergedNewColors)
            if verbose > 3:
                print("   Changing original colors:", flush=True)
            rank = 0
            for color in range(len(SortedRawModuleColors)):
                if SortedRawModuleColors[color] != unassdColor:
                    rank = rank + 1
                    if verbose > 3:
                        print("      ", SortedRawModuleColors[color], "to ", colorSeq[rank], flush=True)
                    MergedNewColors[MergedColors == SortedRawModuleColors[color]] = colorSeq[rank]
            if isinstance(pd.Categorical, MergedColors):
                MergedNewColors = pd.Categorical(MergedNewColors)
        else:
            MergedNewColors = MergedColors
        # MergedNewColors = MergedNewColors[, drop = TRUE]
        if getNewMEs:
            if nNewMods < nOldMods or relabel or getNewUnassdME:
                if verbose > 0:
                    print("  Calculating new MEs...", flush=True)
                NewMEs = multiSetMEs(exprData, colors=None, universalColors=MergedNewColors, useSets=useSets,
                                     impute=impute, subHubs=True, trapErrors=False,
                                     excludeGrey=not getNewUnassdME, grey=unassdColor, verbose=verbose - 1)
                newMEs = consensusOrderMEs(NewMEs, useAbs=useAbs, useSets=useSets, greyLast=True, greyName=greyName)
                ConsDiss = consensusMEDissimilarity(newMEs, equalizeQuantiles=equalizeQuantiles,
                                                    quantileSummary=quantileSummary,
                                                    consensusQuantile=consensusQuantile,
                                                    useAbs=useAbs, useSets=useSets, greyName=greyName)
                if len(ConsDiss) > 1:
                    Tree = hclust(pdist(ConsDiss), method="average")
                else:
                    Tree = None
            else:
                newMEs = MEs
        else:
            newMEs = None
        if MEsInSingleFrame:
            newMEs = newMEs[0]  # $data
            oldMEs = oldMEs[0]  # $data
    except:
        if not trapErrors:
            sys.exit("Error!")
        if verbose > 0:
            print("Warning: merging of modules failed", flush=True)
            print(" --> returning unmerged modules and *no* eigengenes.", flush=True)
        print(f"{bcolors.WARNING}mergeCloseModules: merging of modules failed --> returning unmerged modules and *no* "
              f"eigengenes.\n{bcolors.ENDC}")
        return {"colors": origColors, "allOK": False}
    else:
        return {"colors": MergedNewColors, "dendro": Tree, "oldDendro": oldTree, "cutHeight": cutHeight,
                "oldMEs": oldMEs['data'], "newMEs": newMEs['data'], "allOK": True}


# Take in pearsons r and n (number of experiments) to calculate the t-stat and p value (student's t distribution)
def get_pval(r, n):
    # calculate t-stat, n-2 degrees of freedom
    tstat = r * np.sqrt((n - 2) / (1 - r * r))
    # find p-value for the double-sided test. Students t, n-2 degrees of freedom
    pval = stats.t.sf(np.abs(tstat), n - 2) * 2
    return tstat, pval


def createNetwork(datExpr, networkType="unsigned"):
    intType = networkTypes.index(networkType)
    if intType is None:
        sys.exit(("Unrecognized 'networkType'. Recognized values are", str(networkTypes)))

    genes = datExpr.index.values
    samples = datExpr.columns.values
    datExprCor = np.asmatrix(datExpr.iloc[:, :].corr(method='pearson'))

    if intType == 1:  # signed
        datExprCor = 0.5 + 0.5 * datExprCor

    # Crates graph using the data of the correlation matrix
    G = nx.from_numpy_matrix(datExprCor)

    # relabels the nodes to match the  stocks names
    G = nx.relabel_nodes(G, lambda x: genes[x])

    # shows the edges with their corresponding weights
    G.edges(data=True)

    return datExprCor
