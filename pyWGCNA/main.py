import math
import numpy as np
import pandas as pd
import psutil
import scipy
import scipy.stats as stats
import statistics
import sys
import warnings
from scipy.cluster.hierarchy import linkage, cut_tree
import networkx as nx
import matplotlib.pyplot as plt
from statsmodels.formula.api import ols

# public values
networkTypes = ["unsigned", "signed"]
adjacencyTypes = ["unsigned", "signed"]


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
def goodGenesFun(datExpr, weights=None, useSamples=None, useGenes=None, minFraction=1/2,
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


# Determine cluster under the line # TODO
def cutree(sampleTree, cutHeight=50000):
    cutTree = cut_tree(sampleTree, height=cutHeight)

    return cutTree


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
        maxAlloc = psutil.virtual_memory()[2]
    else:
        maxAlloc = maxMemoryAllocation / 8

    print(maxAlloc)
    maxAlloc = maxAlloc / overheadFactor

    print(maxAlloc)

    if rectangularBlocks:
        blockSz = math.floor(maxAlloc / matrixSize)
    else:
        blockSz = math.floor(math.sqrt(maxAlloc))

    print(matrixSize, blockSz)
    return min(matrixSize, blockSz)


# Calculation of fitting statistics for evaluating scale free topology fit.
def scaleFreeFitIndex(k, nBreaks=10, removeFirst=False):
    discretized_k = pd.cut(k, nBreaks)
    dk = discretized_k.mean()  # tapply(k, discretized_k, mean)
    p_dk = len(discretized_k) / len(k)  # as.vector(tapply(k, discretized.k, length)/length(k))
    breaks1 = range(min(k), max(k), nBreaks + 1)
    hist1 = plt.hist(k, breaks=breaks1, plot=False, right=True)
    dk2 = hist1['mids']
    if dk.isna():
        dk = dk2
    if dk == 0:
        dk = dk2
    if p_dk.isna():
        p_dk = 0
    log_dk = math.log10(dk)
    if removeFirst:
        p_dk = p_dk[1:]
        log_dk = log_dk[1:]
    log_p_dk = math.log10(p_dk + 1e-09)
    df = pd.DataFrame(np.array([log_p_dk, log_dk]), columns=['log_p_dk', 'log_dk'])
    model = ols(formula='log_p_dk ~ log_dk', data=df).fit()
    # lm2 = lm(log.p.dk ~ log.dk + I(10^log.dk))
    datout = pd.DataFrame({'Rsquared.SFT': [model.summary().rsquared],
                           'slope.SFT': [model.summary().coefficients],
                           'truncatedExponentialAdjRsquared': [model.summary().adj.r.squared]})

    return datout


# Call the network topology analysis function
def pickSoftThreshold(data, dataIsExpr=True, weights=None, RsquaredCut=0.85,
                      powerVector=np.concatenate((np.arange(1, 10, 1), np.arange(12, 20, 2))),
                      removeFirst=False, nBreaks=10, blockSize=None, corOptions=None,
                      networkType="unsigned", moreNetworkConcepts=False, gcInterval=None, verbose=0):
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
        blockSize = calBlockSize(nGenes, rectangularBlocks=True, maxMemoryAllocation=2 ^ 30)
        if verbose > 0:
            print("pickSoftThreshold: will use block size ", blockSize, ".", flush=True)

    if gcInterval is None or len(gcInterval) == 0:
        gcInterval = 4 * blockSize

    colname1 = ["Power", "SFT.R.sq", "slope", "truncated R.sq", "mean(k)", "median(k)", "max(k)"]

    if moreNetworkConcepts:
        colname1 = colname1.append(["Density", "Centralization", "Heterogeneity"])

    datout = pd.DataFrame(np.full((len(powerVector), len(colname1)), 666), columns=colname1)
    datout.loc[:, 1] = powerVector
    if verbose > 0:
        print("pickSoftThreshold: calculating connectivity for given powers...\n")
    if verbose == 1:
        pind = ""
    else:
        print("\n")

    datk = np.zeros((nGenes, len(powerVector)))
    nPowers = len(powerVector)
    startG = 1
    lastGC = 0
    #print(data.shape[1])
    #print(data.columns.values)
    if corOptions is None:
        corOptions = pd.DataFrame({'data': data}, index=[range(data.shape[1])])
    else:
        corOptions['x'] = data

    if weights is not None:
        if not dataIsExpr:
            sys.exit("Weights can only be used when 'data' represents expression data ('dataIsExpr' must be TRUE).")
        if data.shape != weights.shape:
            sys.exit("When 'weights' are given, dimensions of 'data' and 'weights' must be the same.")
        corOptions[:, 2] = weights

    while startG <= nGenes:
        endG = min(startG + blockSize - 1, nGenes)
        if verbose > 1:
            print("\n  ..working on genes", startG, "through", endG, "of", nGenes, flush=True)
        nBlockGenes = endG - startG + 1

        useGenes = range(startG, endG)
        nGenes1 = len(useGenes)
        if dataIsExpr:
            corOptions['y'] = data[:, useGenes]
            if weights is not None:
                corOptions['weights.y'] = weights[:, useGenes]
            corx = corOptions.corr(method='pearson')
            if intType == 1:
                corx = abs(corx)
            elif intType == 2:
                corx = (1 + corx) / 2
            elif intType == 3:
                corx[corx < 0] = 0
            if sum(corx.isna()) != 0:
                msg = "Some correlations are NA in block" + startG + ":" + endG + "."
                warnings.warn(msg)
        else:
            corx = data[:, useGenes]

        ind = pd.concat([pd.Series(useGenes), pd.Series(range(len(useGenes)))], axis=1)
        corx[ind] = 1
        datk.local = np.empty((nGenes1, nPowers))
        corxPrev = np.ones(corx.shape)
        powerVector1 = range(0, powerVector.head(n=-1))
        powerSteps = powerVector - powerVector1
        uniquePowerSteps = np.unique(powerSteps)

        def fun(p):
            return corx ^ p

        corxPowers = uniquePowerSteps.applymap(fun)
        corxPowers.columns = uniquePowerSteps
        for j in range(nPowers):
            corxCur = corxPrev * corxPowers[[chr(powerSteps[j])]]
            datk.local[:, j] = corxCur.sum(axis=0) - 1
            corxPrev = corxCur

        datk[startG:endG, :] = datk.local

        startG = endG + 1
        if 0 < gcInterval < startG - lastGC:
            lastGC = startG

        # if verbose == 1:
        #    pind = updateProgInd(endG / nGenes, pind)

    if verbose == 1:
        print("", flush=True)

    for i in range(len(powerVector)):
        khelp = datk[:, i]

        SFT1 = scaleFreeFitIndex(k=khelp, nBreaks=nBreaks, removeFirst=removeFirst)
        datout[i, 2] = SFT1['Rsquared.SFT']
        datout[i, 3] = SFT1['slope.SFT']
        datout[i, 4] = SFT1['truncatedExponentialAdjRsquared']
        datout[i, 5] = statistics.mean(khelp)
        datout[i, 6] = statistics.median(khelp)
        datout[i, 7] = max(khelp)

        if moreNetworkConcepts:
            Density = sum(khelp) / (nGenes * (nGenes - 1))
            datout[i, 8] = Density
            Centralization = nGenes * (max(khelp) - statistics.mean(khelp)) / ((nGenes - 1) * (nGenes - 2))
            datout[i, 9] = Centralization
            Heterogeneity = math.sqrt(nGenes * sum(khelp ^ 2) / sum(khelp) ^ 2 - 1)
            datout[i, 10] = Heterogeneity

    print(round(data.frame(datout), 3))
    ind1 = datout[:, 2] > RsquaredCut
    indcut = None
    if sum(ind1) > 0:
        indcut = min(range(len(ind1))[ind1])
    powerEstimate = powerVector[indcut][[1]]

    return powerEstimate, pd.DataFrame(datout)


def adjacency(datExpr, selectCols=None, type="unsigned", power=6, corOptions=pd.DataFrame(), weights=None,
              distFnc="dist", distOptions="method = 'euclidean'", weightArgNames=["weights.x", "weights.y"]):
    intType = adjacencyTypes.index(type)
    if intType is None:
        msg = "Unrecognized 'type'. Recognized values are", adjacencyTypes
        sys.exit(msg)
    checkAndScaleWeights(weights, datExpr, scaleByMax=False)
    if len(weights) > 0:
        if selectCols.isnull():
            if isinstance(corOptions, pd.DataFrame):
                weightOpt = weights.x = weights
                weightOpt.index = weightArgNames[1]
            else:
                weightOpt = weightArgNames[1] + " = weights"
        else:
            if isinstance(corOptions, pd.DataFrame):
                weightOpt = weights.x = weights, weights.y = weights[:, selectCols]
                weightOpt.index = weightArgNames[1:2]
            else:
                weightOpt = weightArgNames[1] + " = weights, " + weightArgNames[2] + " = weights[, selectCols]"
    else:
        if isinstance(corOptions, pd.DataFrame):
            weightOpt = pd.DataFrame()
        else:
            weightOpt = ""

    if intType < 4:
        if selectCols.isnull():
            if isinstance(corOptions, pd.DataFrame):
                cor_mat = scipy.stats.pearsonr(datExpr, weightOpt, corOptions)
            else:
                corExpr = parse(text=corFnc + "(datExpr " + prepComma(weightOpt) + prepComma(corOptions) + ")")
                cor_mat = eval(corExpr)
        else:
            if isinstance(corOptions, pd.DataFrame):
                cor_mat = scipy.stats.pearsonr(x=datExpr, y=datExpr[:, selectCols])  # , weightOpt, corOptions)
            else:
                corExpr = parse(text=corFnc + "(datExpr, datExpr[, selectCols] " + prepComma(weightOpt) + prepComma(
                    corOptions) + ")")
                cor_mat = eval(corExpr)
    else:
        if not isinstance(selectCols, pd.DataFrame):
            sys.stop("The argument 'selectCols' cannot be used for distance adjacency.")
        if isinstance(distOptions, pd.DataFrame):
            d = scipy.spatial.distance_matrix(datExpr.transpose(), distOptions)
        else:
            corExpr = parse(text=distFnc + "(t(datExpr) " + prepComma(distOptions) + ")")
            d = eval(corExpr)
        if any(d < 0):
            warnings.WarningMessage("Function WGCNA::adjacency: Distance function returned (some) negative values.")
            cor_mat = 1 - ((d / max(d)) ^ 2)

    if intType == 1:
        cor_mat = abs(cor_mat)
    elif intType == 2:
        cor_mat = (1 + cor_mat) / 2
    elif intType == 3:
        cor_mat[cor_mat < 0] = 0

    return cor_mat ^ power


# Take in pearsons r and n (number of experiments) to calculate the t-stat and p value (student's t distribution)
def get_pval(r, n):
    # calculate t-stat, n-2 degrees of freedom
    tstat = r * np.sqrt((n - 2) / (1 - r * r))
    # find p-value for the double-sided test. Students t, n-2 degrees of freedom
    pval = stats.t.sf(np.abs(tstat), n - 2) * 2
    return tstat, pval
