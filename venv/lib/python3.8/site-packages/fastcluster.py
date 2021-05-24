# -*- coding: utf-8 -*-
__doc__ = """Fast hierarchical clustering routines for R and Python

Copyright:
Until package version 1.1.23: © 2011 Daniel Müllner <http://danifold.net>
All changes from version 1.1.24 on: © Google Inc. <http://google.com>

This module provides fast hierarchical clustering routines. The "linkage"
method is designed to provide a replacement for the “linkage” function and
its siblings in the scipy.cluster.hierarchy module. You may use the methods
in this module with the same syntax as the corresponding SciPy functions but
with the benefit of much faster performance.

The method "linkage_vector" performs clustering of vector data with memory-
saving algorithms.

Refer to the User's manual "fastcluster.pdf" for comprehensive details. It
is located in the directory inst/doc/ in the source distribution and may
also be obtained at <http://danifold.net/fastcluster.html>.
"""

__all__ = ['single', 'complete', 'average', 'weighted', 'ward', 'centroid', 'median', 'linkage', 'linkage_vector']
__version_info__ = ('1', '1', '28')
__version__ = '.'.join(__version_info__)

from numpy import double, empty, array, ndarray, var, cov, dot, expand_dims, \
    ceil, sqrt
from numpy.linalg import inv
try:
    from scipy.spatial.distance import pdist
except ImportError:
    def pdist(*args, **kwargs):
        raise ImportError('The fastcluster.linkage function cannot process '
                          'vector data since the function '
                          'scipy.spatial.distance.pdist could not be  '
                          'imported.')
from _fastcluster import linkage_wrap, linkage_vector_wrap

def single(D):
    '''Single linkage clustering (alias). See the help on the “linkage”
function for further information.'''
    return linkage(D, method='single')

def complete(D):
    '''Complete linkage clustering (alias). See the help on the “linkage”
function for further information.'''
    return linkage(D, method='complete')

def average(D):
    '''Hierarchical clustering with the “average” distance update formula
(alias). See the help on the “linkage” function for further information.'''
    return linkage(D, method='average')

def weighted(D):
    '''Hierarchical clustering with the “weighted” distance update formula
(alias). See the help on the “linkage” function for further information.'''
    return linkage(D, method='weighted')

def ward(D):
    '''Hierarchical clustering with the “Ward” distance update formula
(alias). See the help on the “linkage” function for further information.'''
    return linkage(D, method='ward')

def centroid(D):
    '''Hierarchical clustering with the “centroid” distance update formula
(alias). See the help on the “linkage” function for further information.'''
    return linkage(D, method='centroid')

def median(D):
    '''Hierarchical clustering with the “median” distance update formula
(alias). See the help on the “linkage” function for further information.'''
    return linkage(D, method='median')

# This dictionary must agree with the enum method_codes in fastcluster.cpp.
mthidx = {'single'   : 0,
          'complete' : 1,
          'average'  : 2,
          'weighted' : 3,
          'ward'     : 4,
          'centroid' : 5,
          'median'   : 6 }

def linkage(X, method='single', metric='euclidean', preserve_input=True):
    r'''Hierarchical, agglomerative clustering on a dissimilarity matrix or on
Euclidean data.

Apart from the argument 'preserve_input', the method has the same input
parameters and output format as the functions of the same name in the
module scipy.cluster.hierarchy.

The argument X is preferably a NumPy array with floating point entries
(X.dtype==numpy.double). Any other data format will be converted before
it is processed.

If X is a one-dimensional array, it is considered a condensed matrix of
pairwise dissimilarities in the format which is returned by
scipy.spatial.distance.pdist. It contains the flattened, upper-
triangular part of a pairwise dissimilarity matrix. That is, if there
are N data points and the matrix d contains the dissimilarity between
the i-th and j-th observation at position d(i,j), the vector X has
length N(N-1)/2 and is ordered as follows:

  [ d(0,1), d(0,2), ..., d(0,n-1), d(1,2), ..., d(1,n-1), ...,
    d(n-2,n-1) ]

The 'metric' argument is ignored in case of dissimilarity input.

The optional argument 'preserve_input' specifies whether the method
makes a working copy of the dissimilarity vector or writes temporary
data into the existing array. If the dissimilarities are generated for
the clustering step only and are not needed afterward, approximately
half the memory can be saved by specifying 'preserve_input=False'. Note
that the input array X contains unspecified values after this procedure.
It is therefore safer to write

  linkage(X, method="...", preserve_input=False)
  del X

to make sure that the matrix X is not accessed accidentally after it has
been used as scratch memory. (The single linkage algorithm does not
write to the distance matrix or its copy anyway, so the 'preserve_input'
flag has no effect in this case.)

If X contains vector data, it must be a two-dimensional array with N
observations in D dimensions as an (N×D) array. The preserve_input
argument is ignored in this case. The specified metric is used to
generate pairwise distances from the input. The following two function
calls yield the same output:

  linkage(pdist(X, metric), method="...", preserve_input=False)
  linkage(X, metric=metric, method="...")

The general scheme of the agglomerative clustering procedure is as
follows:

  1. Start with N singleton clusters (nodes) labeled 0,...,N−1, which
     represent the input points.
  2. Find a pair of nodes with minimal distance among all pairwise
     distances.
  3. Join the two nodes into a new node and remove the two old nodes.
     The new nodes are labeled consecutively N, N+1, ...
  4. The distances from the new node to all other nodes is determined by
     the method parameter (see below).
  5. Repeat N−1 times from step 2, until there is one big node, which
     contains all original input points.

The output of linkage is stepwise dendrogram, which is represented as an
(N−1)×4 NumPy array with floating point entries (dtype=numpy.double).
The first two columns contain the node indices which are joined in each
step. The input nodes are labeled 0,...,N−1, and the newly generated
nodes have the labels N,...,2N−2. The third column contains the distance
between the two nodes at each step, ie. the current minimal distance at
the time of the merge. The fourth column counts the number of points
which comprise each new node.

The parameter method specifies which clustering scheme to use. The
clustering scheme determines the distance from a new node to the other
nodes. Denote the dissimilarities by d, the nodes to be joined by I, J,
the new node by K and any other node by L. The symbol |I| denotes the
size of the cluster I.

method='single': d(K,L) = min(d(I,L), d(J,L))

  The distance between two clusters A, B is the closest distance between
  any two points in each cluster:

    d(A,B) = min{ d(a,b) | a∈A, b∈B }

method='complete': d(K,L) = max(d(I,L), d(J,L))

  The distance between two clusters A, B is the maximal distance between
  any two points in each cluster:

    d(A,B) = max{ d(a,b) | a∈A, b∈B }

method='average': d(K,L) = ( |I|·d(I,L) + |J|·d(J,L) ) / (|I|+|J|)

  The distance between two clusters A, B is the average distance between
  the points in the two clusters:

    d(A,B) = (|A|·|B|)^(-1) · \sum { d(a,b) | a∈A, b∈B }

method='weighted': d(K,L) = (d(I,L)+d(J,L))/2

  There is no global description for the distance between clusters since
  the distance depends on the order of the merging steps.

The following three methods are intended for Euclidean data only, ie.
when X contains the pairwise (non-squared!) distances between vectors in
Euclidean space. The algorithm will work on any input, however, and it
is up to the user to make sure that applying the methods makes sense.

method='centroid': d(K,L) = ( (|I|·d(I,L) + |J|·d(J,L)) / (|I|+|J|)
                              − |I|·|J|·d(I,J)/(|I|+|J|)^2 )^(1/2)

  There is a geometric interpretation: d(A,B) is the distance between
  the centroids (ie. barycenters) of the clusters in Euclidean space:

    d(A,B) = ‖c_A−c_B∥,

  where c_A denotes the centroid of the points in cluster A.

method='median': d(K,L) = ( d(I,L)/2 + d(J,L)/2 − d(I,J)/4 )^(1/2)

  Define the midpoint w_K of a cluster K iteratively as w_K=k if K={k}
  is a singleton and as the midpoint (w_I+w_J)/2 if K is formed by
  joining I and J. Then we have

    d(A,B) = ∥w_A−w_B∥

  in Euclidean space for all nodes A,B. Notice however that this
  distance depends on the order of the merging steps.

method='ward': d(K,L) = ( ((|I|+|L)d(I,L) + (|J|+|L|)d(J,L) − |L|d(I,J))
                          / (|I|+|J|+|L|) )^(1/2)

  The global cluster dissimilarity can be expressed as

    d(A,B) = ( 2|A|·|B|/(|A|+|B|) )^(1/2) · ‖c_A−c_B∥,

  where c_A again denotes the centroid of the points in cluster A.

The clustering algorithm handles infinite values correctly, as long as the
chosen distance update formula makes sense. If a NaN value occurs, either
in the original dissimilarities or as an updated dissimilarity, an error is
raised.

The linkage method does not treat NumPy's masked arrays as special
and simply ignores the mask.'''
    X = array(X, copy=False, subok=True)
    if X.ndim==1:
        if method=='single':
            preserve_input = False
        X = array(X, dtype=double, copy=preserve_input, order='C', subok=True)
        NN = len(X)
        N = int(ceil(sqrt(NN*2)))
        if (N*(N-1)//2) != NN:
            raise ValueError(r'The length of the condensed distance matrix '
                             r'must be (k \choose 2) for k data points!')
    else:
        assert X.ndim==2
        N = len(X)
        X = pdist(X, metric=metric)
        X = array(X, dtype=double, copy=False, order='C', subok=True)
    Z = empty((N-1,4))
    if N > 1:
        linkage_wrap(N, X, Z, mthidx[method])
    return Z

# This dictionary must agree with the enum metric_codes in fastcluster_python.cpp.
mtridx = {'euclidean'      :  0,
          'minkowski'      :  1,
          'cityblock'      :  2,
          'seuclidean'     :  3,
          'sqeuclidean'    :  4,
          'cosine'         :  5,
          'hamming'        :  6,
          'jaccard'        :  7,
          'chebychev'      :  8,
          'canberra'       :  9,
          'braycurtis'     : 10,
          'mahalanobis'    : 11,
          'yule'           : 12,
          'matching'       : 13,
          'sokalmichener'  : 13, # an alias for 'matching'
          'dice'           : 14,
          'rogerstanimoto' : 15,
          'russellrao'     : 16,
          'sokalsneath'    : 17,
          'kulsinski'      : 18,
          'USER'           : 19,
          }

booleanmetrics = ('yule', 'matching', 'dice', 'kulsinski', 'rogerstanimoto',
                  'sokalmichener', 'russellrao', 'sokalsneath', 'kulsinski')

def linkage_vector(X, method='single', metric='euclidean', extraarg=None):
    r'''Hierarchical (agglomerative) clustering on Euclidean data.

Compared to the 'linkage' method, 'linkage_vector' uses a memory-saving
algorithm. While the linkage method requires Θ(N^2) memory for
clustering of N points, this method needs Θ(ND) for N points in R^D,
which is usually much smaller.

The argument X has the same format as before, when X describes vector
data, ie. it is an (N×D) array. Also the output array has the same
format. The parameter method must be one of 'single', 'centroid',
'median', 'ward', ie. only for these methods there exist memory-saving
algorithms currently. If 'method', is one of 'centroid', 'median',
'ward', the 'metric' must be 'euclidean'.

For single linkage clustering, any dissimilarity function may be chosen.
Basically, every metric which is implemented in the method
scipy.spatial.distance.pdist is reimplemented here. However, the metrics
differ in some instances since a number of mistakes and typos (both in
the code and in the documentation) were corrected in the fastcluster
package.

Therefore, the available metrics with their definitions are listed below
as a reference. The symbols u and v mostly denote vectors in R^D with
coordinates u_j and v_j respectively. See below for additional metrics
for Boolean vectors. Unless otherwise stated, the input array X is
converted to a floating point array (X.dtype==numpy.double) if it does
not have already the required data type. Some metrics accept Boolean
input; in this case this is stated explicitly below.

If a NaN value occurs, either in the original dissimilarities or as an
updated dissimilarity, an error is raised. In principle, the clustering
algorithm handles infinite values correctly, but the user is advised to
carefully check the behavior of the metric and distance update formulas
under these circumstances.

The distance formulas combined with the clustering in the
'linkage_vector' method do not have specified behavior if the data X
contains infinite or NaN values. Also, the masks in NumPy’s masked
arrays are simply ignored.

metric='euclidean': Euclidean metric, L_2 norm

    d(u,v) = ∥u−v∥ = ( \sum_j { (u_j−v_j)^2 } )^(1/2)

metric='sqeuclidean': squared Euclidean metric

    d(u,v) = ∥u−v∥^2 = \sum_j { (u_j−v_j)^2 }

metric='seuclidean': standardized Euclidean metric

    d(u,v) = ( \sum_j { (u_j−v_j)^2 / V_j } )^(1/2)

  The vector V=(V_0,...,V_{D−1}) is given as the 'extraarg' argument. If
  no 'extraarg' is given, V_j is by default the unbiased sample variance
  of all observations in the j-th coordinate:

    V_j = Var_i (X(i,j) ) = 1/(N−1) · \sum_i ( X(i,j)^2 − μ(X_j)^2 )

  (Here, μ(X_j) denotes as usual the mean of X(i,j) over all rows i.)

metric='mahalanobis': Mahalanobis distance

    d(u,v) = ( transpose(u−v) V (u−v) )^(1/2)

  Here, V=extraarg, a (D×D)-matrix. If V is not specified, the inverse
  of the covariance matrix numpy.linalg.inv(numpy.cov(X, rowvar=False))
  is used.

metric='cityblock': the Manhattan distance, L_1 norm

    d(u,v) = \sum_j |u_j−v_j|

metric='chebychev': the supremum norm, L_∞ norm

    d(u,v) = max_j { |u_j−v_j| }

metric='minkowski': the L_p norm

    d(u,v) = ( \sum_j |u_j−v_j|^p ) ^(1/p)

  This metric coincides with the cityblock, euclidean and chebychev
  metrics for p=1, p=2 and p=∞ (numpy.inf), respectively. The parameter
  p is given as the 'extraarg' argument.

metric='cosine'

    d(u,v) = 1 − ⟨u,v⟩ / (∥u∥·∥v∥)
           = 1 − (\sum_j u_j·v_j) / ( (\sum u_j^2)(\sum v_j^2) )^(1/2)

metric='correlation': This method first mean-centers the rows of X and
  then applies the 'cosine' distance. Equivalently, the correlation
  distance measures 1 − (Pearson’s correlation coefficient).

    d(u,v) = 1 − ⟨u−μ(u),v−μ(v)⟩ / (∥u−μ(u)∥·∥v−μ(v)∥)

metric='canberra'

    d(u,v) = \sum_j ( |u_j−v_j| / (|u_j|+|v_j|) )

  Summands with u_j=v_j=0 contribute 0 to the sum.

metric='braycurtis'

    d(u,v) = (\sum_j |u_j-v_j|) / (\sum_j |u_j+v_j|)

metric=(user function): The parameter metric may also be a function
  which accepts two NumPy floating point vectors and returns a number.
  Eg. the Euclidean distance could be emulated with

    fn = lambda u, v: numpy.sqrt(((u-v)*(u-v)).sum())
    linkage_vector(X, method='single', metric=fn)

  This method, however, is much slower than the build-in function.

metric='hamming': The Hamming distance accepts a Boolean array
  (X.dtype==bool) for efficient storage. Any other data type is
  converted to numpy.double.

    d(u,v) = |{j | u_j≠v_j }|

metric='jaccard': The Jaccard distance accepts a Boolean array
  (X.dtype==bool) for efficient storage. Any other data type is
  converted to numpy.double.

    d(u,v) = |{j | u_j≠v_j }| / |{j | u_j≠0 or v_j≠0 }|
    d(0,0) = 0

  Python represents True by 1 and False by 0. In the Boolean case, the
  Jaccard distance is therefore:

    d(u,v) = |{j | u_j≠v_j }| / |{j | u_j  ∨ v_j }|

The following metrics are designed for Boolean vectors. The input array
is converted to the 'bool' data type if it is not Boolean already. Use
the following abbreviations to count the number of True/False
combinations:

  a = |{j | u_j ∧ v_j }|
  b = |{j | u_j ∧ (¬v_j) }|
  c = |{j | (¬u_j) ∧ v_j }|
  d = |{j | (¬u_j) ∧ (¬v_j) }|

Recall that D denotes the number of dimensions, hence D=a+b+c+d.

metric='yule'

    d(u,v) = 2bc / (ad+bc)

metric='dice':

    d(u,v) = (b+c) / (2a+b+c)
    d(0,0) = 0

metric='rogerstanimoto':

    d(u,v) = 2(b+c) / (b+c+D)

metric='russellrao':

    d(u,v) = (b+c+d) / D

metric='sokalsneath':

    d(u,v) = 2(b+c)/ ( a+2(b+c))
    d(0,0) = 0

metric='kulsinski'

    d(u,v) = (b/(a+b) + c/(a+c)) / 2

metric='matching':

    d(u,v) = (b+c)/D

  Notice that when given a Boolean array, the 'matching' and 'hamming'
  distance are the same. The 'matching' distance formula, however,
  converts every input to Boolean first. Hence, the vectors (0,1) and
  (0,2) have zero 'matching' distance since they are both converted to
  (False, True) but the Hamming distance is 0.5.

metric='sokalmichener' is an alias for 'matching'.'''
    if method=='single':
        assert metric!='USER'
        if metric in ('hamming', 'jaccard'):
            X = array(X, copy=False, subok=True)
            dtype = bool if X.dtype==bool else double
        else:
            dtype = bool if metric in booleanmetrics else double
        X = array(X, dtype=dtype, copy=False, order='C', subok=True)
    else:
        assert metric=='euclidean'
        X = array(X, dtype=double, copy=(method=='ward'), order='C', subok=True)
    assert X.ndim==2
    N = len(X)
    Z = empty((N-1,4))

    if metric=='seuclidean':
        if extraarg is None:
            extraarg = var(X, axis=0, ddof=1)
    elif metric=='mahalanobis':
        if extraarg is None:
            extraarg = inv(cov(X, rowvar=False))
        # instead of the inverse covariance matrix, pass the matrix product
        # with the data matrix!
        extraarg = array(dot(X,extraarg),dtype=double, copy=False, order='C', subok=True)
    elif metric=='correlation':
        X = X-expand_dims(X.mean(axis=1),1)
        metric='cosine'
    elif not isinstance(metric, str):
        assert extraarg is None
        metric, extraarg = 'USER', metric
    elif metric!='minkowski':
        assert extraarg is None
    if N > 1:
        linkage_vector_wrap(X, Z, mthidx[method], mtridx[metric], extraarg)
    return Z
