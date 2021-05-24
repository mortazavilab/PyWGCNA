import math, numpy as np, pandas as pd
import sys, psutil, warnings
import fastcluster


#public values
import numpy

networkTypes = ["unsigned", "signed", "signed hybrid"]

def checkAndScaleWeights(weights, expr, scaleByMax=True, verbose=1):
    if len(weights) == 0:
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
def goodGenes(datExpr, weights=None, useSamples=None, useGenes=None, minFraction=1/2,
              minNSamples=4, minNGenes=4, tol=None, minRelativeWeight=0.1, verbose=1):
    datExpr = np.asmatrix(datExpr)
    if not (all(ele.isdigit() for ele in datExpr)):
        sys.exit("datExpr must contain numeric data.")

    weights =checkAndScaleWeights(weights, datExpr, scaleByMax=True)
    if tol is None:
        tol = 1e-10 * max(abs(datExpr))
    if useGenes is None:
        useGenes = np.repeat(True, datExpr.shape[1])
    if useSamples is None:
        useSamples = np.repeat(True, datExpr.shape[0])

    if len(useGenes) != datExpr.shape[1]:
        sys.exit("Length of nGenes is not compatible with number of columns in datExpr.")
    if len(useSamples) != datExpr.shape[0]:
        sys.exit("Length of nSamples is not compatible with number of rows in datExpr.")

    nSamples = sum(useSamples)
    nGenes = sum(useGenes)
    if len(weights) == 0:
        nPresent = datExpr[datExpr[useSamples, useGenes].notna()].sum(axis = 0)
    else :
        nPresent = datExpr[datExpr[useSamples, useGenes].notna() and
                           weights[useSamples, useGenes] > minRelativeWeight].sum(axis = 0)

    gg = useGenes
    gg[useGenes][nPresent < minNSamples] = False

    if len(weights) == 0:
        var = datExpr[useSamples, gg].var(axis = 0)
    else:
        ## need to be fix
        #TODO:colWeightedVars
        var = np.var(datExpr, w = weights)

    var[var.isna()] = 0
    nNAsGenes = datExpr[useSamples, gg].isna().sum(axis=0)
    gg[gg] = (nNAsGenes < (1 - minFraction) * nSamples & var > tol ^ 2 and
              (nSamples - nNAsGenes >= minNSamples))

    if sum(gg) < minNGenes:
        sys.exit("Too few genes with valid expression levels in the required number of samples.")
    if verbose > 0 and nGenes - sum(gg) > 0:
        print("\n\n  ..Excluding", nGenes - sum(gg),
              "genes from the calculation due to too many missing samples or zero variance.\n\n")

    return gg


# Check that all genes and samples have sufficiently low numbers of missing values.
def goodSamplesGenes(datExpr, weights=None, minFraction=1/2, minNSamples=4, minNGenes=4,
                      tol=None, minRelativeWeight=0.1, verbose=1, indent=0):
    goodGenes = None
    goodSamples = None
    nBadGenes = 0
    nBadSamples = 0
    changed = True
    iter = 1
    if (verbose > 0):
        print("Flagging genes and samples with too many missing values...\n\n")
    while (changed):
        if (verbose > 0):
            print(" ..step", iter, "\n")
        goodGenes = goodGenes(datExpr, weights, goodSamples, goodGenes, minFraction = minFraction,
                              minNSamples = minNSamples, minNGenes = minNGenes,
                              minRelativeWeight = minRelativeWeight, tol = tol, verbose = verbose - 1,
                              indent = indent + 1)
        goodSamples = goodSamples(datExpr, weights, goodSamples, goodGenes, minFraction = minFraction,
                                  minNSamples = minNSamples, minNGenes = minNGenes,
                                  minRelativeWeight = minRelativeWeight, verbose = verbose - 1,
                                  indent = indent + 1)
        changed = (~goodGenes.sum() > nBadGenes) or (~goodSamples.sum() > nBadSamples)
        nBadGenes = ~goodGenes.sum()
        nBadSamples = ~goodSamples.sum()
        iter = iter + 1

    allOK = np.concatenate(nBadGenes, nBadSamples).sum() == 0

    return goodGenes, goodSamples, allOK


def hclust(d, method = "complete", members = None):
    if method == "ward":
        print("\nThe \"ward\" method has been renamed to \"ward.D\"; note new \"ward.D2\"\n")
        method = "ward.D"

    METHODS = ["single", "complete", "average", "mcquitty", "ward.D", "centroid", "median", "ward.D2"]

    if method not in METHODS:
        sys.exit("Invalid clustering method.")

    if method == -1:
        sys.exit("Ambiguous clustering method.")

    dendrogram = fastcluster.linkage(d, method=method)

    return dendrogram


# Determine cluster under the line # TODO
def cutree(sampleTree, cutHeight=50000):
    return sampleTree


def checkSimilarity(adjMat, min = 0, max = 1):
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


def calBlockSize(matrixSize, rectangularBlocks = True, maxMemoryAllocation = None, overheadFactor = 3):
    if maxMemoryAllocation is None:
        maxAlloc = psutil.virtual_memory()[2]
    else:
        maxAlloc = maxMemoryAllocation / 8

    maxAlloc = maxAlloc / overheadFactor

    if rectangularBlocks:
        blockSz = math.floor(maxAlloc / matrixSize)
    else:
        blockSz = math.floor(math.sqrt(maxAlloc))

    return min(matrixSize, blockSz)


# Call the network topology analysis function
def pickSoftThreshold(data, dataIsExpr=True, weights=None, RsquaredCut=0.85,
                      powerVector = np.concatenate((np.arange(1, 10, 1), np.arange(12, 20, 2))),
                      removeFirst=False, nBreaks = 10, blockSize=None, corOptions = list(use = "p"),
                      networkType = "unsigned",  moreNetworkConcepts = False, gcInterval = None, verbose = 0):

    powerVector = np.sort(powerVector)
    intType = networkTypes.index(networkType)
    if intType is None:
        sys.exit(("Unrecognized 'networkType'. Recognized values are", str(networkTypes)))

    nGenes = data.shape[1] #col
    if nGenes < 3:
        sys.exit("The input data data contain fewer than 3 rows (nodes).\n"
                 "This would result in a trivial correlation network.")

    if ~dataIsExpr:
        checkSimilarity(data)
    if any(np.diag(data) != 1):
        data = np.where(np.diag(data), 1)

    if blockSize is None:
        blockSize = calBlockSize(nGenes, rectangularBlocks = True, maxMemoryAllocation = 2 ^ 30)

    if verbose > 0:
        print("pickSoftThreshold: will use block size ", blockSize, ".", flush=True)
    if len(gcInterval) == 0:
        gcInterval = 4 * blockSize

    colname1 = ["Power", "SFT.R.sq", "slope", "truncated R.sq", "mean(k)", "median(k)", "max(k)"]

    if moreNetworkConcepts:
        colname1 = [colname1, "Density", "Centralization", "Heterogeneity"]

    datout = pd.DataFrame(np.full((len(powerVector), len(colname1)), 666), columns=colname1)
    datout[:, 1] = powerVector
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
    corOptions$x = data

    if weights is not None:
        if not dataIsExpr:
            sys.exit("Weights can only be used when 'data' represents expression data ('dataIsExpr' must be TRUE).")
        if data.shape != weights.shape:
            sys.exit("When 'weights' are given, dimensions of 'data' and 'weights' must be the same.")
        corOptions$weights.x = weights

    while startG <= nGenes:
        endG = min(startG + blockSize - 1, nGenes)
        if verbose > 1:
            print("\n  ..working on genes", startG, "through", endG, "of", nGenes, flush=True)
        nBlockGenes = endG - startG + 1
        datk[startG:endG, :] = foreach(t=actualThreads, .combine = rbind) % dopar % {useGenes = c(startG:endG)[jobs[[t]]]
        nGenes1 = len(useGenes)
        if dataIsExpr:
            corOptions$y = data[:, useGenes]
            if weights is not None:
                corOptions$weights.y = weights[, useGenes]
            corx = do.call(corFnc, corOptions)
            if intType == 1:
                corx = math.abs(corx)
            elif intType == 2:
                corx = (1 + corx) / 2
            elif intType == 3:
                corx[corx < 0] = 0
            if sum(corx.isna()) != 0:
                warnings.WarningMessage("Some correlations are NA in block", startG, ":", endG, ".")
        else:
            corx = data[:, useGenes]

        ind = cbind(useGenes, 1:len(useGenes))
        corx[ind] = 1
    datk.local = matrix(NA, nGenes1, nPowers)
    corxPrev = matrix(1, nrow=nrow(corx), ncol=ncol(corx))
    powerVector1 < - c(0, head(powerVector, -1))
    powerSteps < - powerVector - powerVector1
    uniquePowerSteps < - unique(powerSteps)
    corxPowers < - lapply(uniquePowerSteps, function(p)
    corx ^ p)
    names(corxPowers) < - uniquePowerSteps
    for (j in 1:nPowers) {
        corxCur < - corxPrev * corxPowers[[as.character(powerSteps[j])]]
        datk.local[, j] = colSums(corxCur, na.rm = TRUE) - 1
        corxPrev < - corxCur

startG = endG + 1
if ((gcInterval > 0) & & (startG - lastGC > gcInterval)) {
gc()
lastGC = startG
}
if (verbose == 1)
pind = updateProgInd(endG / nGenes, pind)
}
if (verbose == 1)
    printFlush("")
for (i in c(1:length(powerVector))) {
khelp = datk[, i]
if (any(khelp < 0))
browser()
SFT1 = scaleFreeFitIndex(k = khelp, nBreaks = nBreaks,
removeFirst = removeFirst)
datout[i, 2] = SFT1$Rsquared.SFT
datout[i, 3] = SFT1$slope.SFT
datout[i, 4] = SFT1$truncatedExponentialAdjRsquared
datout[i, 5] = mean(khelp, na.rm = TRUE)
datout[i, 6] = median(khelp, na.rm = TRUE)
datout[i, 7] = max(khelp, na.rm = TRUE)
if (moreNetworkConcepts) {
Density = sum(khelp) / (nGenes * (nGenes - 1))
datout[i, 8] = Density
Centralization = nGenes * (max(khelp) - mean(khelp)) / ((nGenes -
1) * (nGenes - 2))
datout[i, 9] = Centralization
Heterogeneity = sqrt(nGenes * sum(khelp ^ 2) / sum(khelp) ^ 2 -
1)
datout[i, 10] = Heterogeneity
}
}
print(signif(data.frame(datout), 3))
ind1 = datout[, 2] > RsquaredCut
indcut = NA
indcut =
if (sum(ind1) > 0)
    min(c(1: length(ind1))[ind1])
    else indcut
    powerEstimate = powerVector[indcut][[1]]
    gc()
    list(powerEstimate=powerEstimate, fitIndices=data.frame(datout))
