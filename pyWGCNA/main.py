import math
import numpy as np
import pandas as pd
import scipy
import scipy.stats as stats
import statistics
import sys
import warnings
from scipy.cluster.hierarchy import linkage, cut_tree
import networkx as nx
from statsmodels.formula.api import ols
import resource

# public values
networkTypes = ["unsigned", "signed"]
adjacencyTypes = ["unsigned", "signed"]
TOMTypes = ["NA", "unsigned", "signed"]
TOMDenoms = ["min", "mean"]

# remove runtime warning (divided by zero)
np.seterr(divide='ignore', invalid='ignore')


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
            warnings.WarningMessage("Found non-finite weights. The corresponding data points will be removed.")
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
def cutree(sampleTree, cutHeight=50000):
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
def pickSoftThreshold(data, dataIsExpr=True, weights=None,
                      powerVector=np.concatenate((np.arange(1, 10, 1), np.arange(12, 20, 2))),
                      nBreaks=10, blockSize=None, corOptions=None, networkType="unsigned",
                      moreNetworkConcepts=False, gcInterval=None, verbose=0):
    powerVector = np.sort(powerVector)
    intType = networkTypes.index(networkType)
    if intType is None:
        sys.exit(("Unrecognized 'networkType'. Recognized values are", str(networkTypes)))

    nGenes = data.shape[1]  # col
    if nGenes < 3:
        sys.exit("The input data data contain fewer than 3 rows (nodes).\n"
                 "This would result in a trivial correlation network.")

    if not dataIsExpr:
        checkSimilarity(data)
        if any(np.diag(data) != 1):
            data = np.where(np.diag(data), 1)

    if blockSize is None:
        blockSize = calBlockSize(nGenes, rectangularBlocks=True, maxMemoryAllocation=2 ** 30)
        if verbose > 0:
            print("pickSoftThreshold: will use block size ", blockSize, flush=True)

    if gcInterval is None or len(gcInterval) == 0:
        gcInterval = 4 * blockSize

    colname1 = ["Power", "SFT.R.sq", "slope", "truncated R.sq", "mean(k)", "median(k)", "max(k)"]

    if moreNetworkConcepts:
        colname1 = colname1.append(["Density", "Centralization", "Heterogeneity"])

    datout = pd.DataFrame(np.full((len(powerVector), len(colname1)), 666), columns=colname1, dtype=object)
    datout['Power'] = powerVector
    if verbose > 0:
        print("pickSoftThreshold: calculating connectivity for given powers...", flush=True)
    else:
        print("\n")

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
        corOptions[:, 2] = weights

    while startG <= nGenes:
        endG = min(startG + blockSize, nGenes)
        if verbose > 1:
            print("\t..working on genes", startG, "through", endG, "of", nGenes, flush=True)

        useGenes = list(range(startG, endG))
        nGenes1 = len(useGenes)
        if dataIsExpr:
            corOptions['y'] = [data.iloc[:, useGenes]]
            if weights is not None:
                corOptions['weights.y'] = [weights.iloc[:, useGenes]]
            corx = np.corrcoef(corOptions.x[0], corOptions.y[0], rowvar=False)
            corx = corx[0:corOptions.x[0].shape[1], 0:corOptions.y[0].shape[1]]
            if intType == 0:
                corx = abs(corx)
            elif intType == 1:
                corx = (1 + corx) / 2

            if np.count_nonzero(np.isnan(corx)) != 0:
                msg = "Some correlations are NA in block " + str(startG) + ":" + str(endG) + "."
                warnings.warn(msg)
        else:
            corx = data[:, useGenes]

        corx[useGenes, list(range(len(useGenes)))] = 1
        datk_local = np.empty((nGenes1, nPowers))
        datk_local[:] = np.NaN
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

        startG = endG + 1
        if 0 < gcInterval < startG - lastGC:
            lastGC = startG

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

    return datout


def adjacency(datExpr, selectCols=None, networkType="unsigned", power=6, corOptions=pd.DataFrame(), weights=None,
              weightArgNames=["weights.x", "weights.y"]):
    intType = networkTypes.index(networkType)
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


def TOMsimilarity(adjMat, TOMType="unsigned", TOMDenom="min",
                  suppressNegativeTOM=False, verbose=1, indent=0):
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

    if verbose > 0:
        print("..connectivity..\n")

    tom = TomSimilarityFromAdj(adjMat, TOMDenomC, TOMTypeC)

    if verbose > 0:
        print("..done..\n")

    return tom


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
