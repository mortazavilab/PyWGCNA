# Quick start: initializing, adding data, runing and saving your PyWGCNA

First, make sure to [install PyWGCNA](https://github.com/mortazavilab/PyWGCNA#readme).

We encourage you to look at [Data input, cleaning and pre-processing] tutorial (Data%20format.md).

Then, download the data from [here](https://github.com/mortazavilab/PyWGCNA/raw/fb2cdc6e4fc1d7ec20ee6e9c39930409018c09c0/tutorials/5xFAD_paper.zip) and uncompressed it.

For this tutorial, we will be analyzing the TPM matrix of 5xFAD mouse from [MODEL-AD](https://www.model-ad.org/) portal and [this paper](https://www.nature.com/articles/s41597-021-01054-y) in Scientific Data which has 5548 genes and 193 samples in four different time point.

## Setupping up PyWGCNA object

We use the TPM matrix to create a PyWGCNA object. The object serves as a container that contains both data (like the TPM matrix) and analysis (like clustering, or visualizing results) for a Bulk RNA-seq dataset. For a technical discussion of the PyWGCNA object structure, check out our [GitHub Wiki](https://github.com/mortazavilab/PyWGCNA.wiki.git). For example, the raw TPM matrix is stored in `GeneExp` class.


```python
import PyWGCNA
geneExp = '5xFAD_paper/expressionList_sorted'
pyWGCNA_5xFAD = PyWGCNA.WGCNA(name='5xFAD', species='mouse',
                              geneExpPath=geneExp, 
                              sep='\t', save=True)
pyWGCNA_5xFAD.expressionList.head(5)
```

    /Users/nargesrezaie/miniconda3/lib/python3.8/site-packages/anndata/_core/anndata.py:120: ImplicitModificationWarning: Transforming to str index.
      warnings.warn("Transforming to str index.", ImplicitModificationWarning)



    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    /var/folders/26/90_v4ll95wjdp9l6g2lxlm8r0000gn/T/ipykernel_39001/2943457033.py in <module>
          1 import PyWGCNA
          2 geneExp = '5xFAD_paper/expressionList_sorted'
    ----> 3 pyWGCNA_5xFAD = PyWGCNA.WGCNA(name='5xFAD', species='mouse',
          4                               geneExpPath=geneExp,
          5                               sep='\t', save=True)


    ~/miniconda3/lib/python3.8/site-packages/PyWGCNA/wgcna.py in __init__(self, name, TPMcutoff, powers, RsquaredCut, MeanCut, networkType, TOMType, minModuleSize, naColor, cut, MEDissThres, species, level, anndata, geneExp, geneExpPath, sep, save, outputPath)
        177                  save=False, outputPath=None):
        178 
    --> 179         super().__init__(species=species, level=level, anndata=anndata, geneExp=geneExp,
        180                          geneExpPath=geneExpPath, sep=sep)
        181 


    ~/miniconda3/lib/python3.8/site-packages/PyWGCNA/geneExp.py in __init__(self, species, level, anndata, geneExp, geneExpPath, sep)
         59         expressionList = expressionList.drop([expressionList.columns[0]], axis=1)
         60 
    ---> 61         self.geneExpr = ad.AnnData(X=expressionList, obs=geneInfo, var=sampleInfo)
         62 
         63     def updateGeneInfo(self, geneInfo=None, path=None, sep=' ', order=True):


    ~/miniconda3/lib/python3.8/site-packages/anndata/_core/anndata.py in __init__(self, X, obs, var, uns, obsm, varm, layers, raw, dtype, shape, filename, filemode, asview, obsp, varp, oidx, vidx)
        306             self._init_as_view(X, oidx, vidx)
        307         else:
    --> 308             self._init_as_actual(
        309                 X=X,
        310                 obs=obs,


    ~/miniconda3/lib/python3.8/site-packages/anndata/_core/anndata.py in _init_as_actual(self, X, obs, var, uns, obsm, varm, varp, obsp, raw, layers, dtype, shape, filename, filemode)
        502                 attr.index = idx
        503             elif not idx.equals(attr.index):
    --> 504                 raise ValueError(f"Index of {attr_name} must match {x_name} of X.")
        505 
        506         # unstructured annotations


    ValueError: Index of obs must match index of X.


## Pre-processing workflow

PyWGCNA allows you to easily preproces the data including removing genes with too many missing values or really low expressed across samples(in default we suggest to remove genes without any expression more than 1 TPM) and also removing samples with too many missing values or not matched with. keep in your mind you can change criteria of removing outlier genes or sample by changing `TPMcutoff` and `cut`


```python
pyWGCNA_5xFAD.preprocess()
```

    [1m[94mPre-processing...[0m
    	Detecting genes and samples with too many missing values...
    	Done pre-processing..
    



    
![png](Quick_Start_files/Quick_Start_3_1.png)
    


## Construction of the gene network and identification of modules

PyWGCNA compress all the steps of network construction and module detection in one function called `findModules` including:
1. Choosing the soft-thresholding power: analysis of network topology
2. Co-expression similarity and adjacency
3. Topological Overlap Matrix (TOM)
4. Clustering using TOM
5. Merging of modules whose expression profiles are very similar


```python
pyWGCNA_5xFAD.findModules()
```

    [1m[94mRun WGCNA...[0m
    [96mpickSoftThreshold: calculating connectivity for given powers...[0m
    will use block size  1876
        Power  SFT.R.sq     slope truncated R.sq      mean(k)    median(k)  \
    0       1  0.368857 -0.481613       0.701585  2444.750756  2260.416614   
    1       2    0.7253  -0.99165       0.886361   840.665489   673.081241   
    2       3  0.791986 -1.194264       0.946969   385.685335   258.451265   
    3       4  0.835392   -1.3419       0.968446   207.404152   113.456087   
    4       5  0.853842 -1.472183       0.973346   123.232581    54.784481   
    5       6  0.870673 -1.553348       0.979584    78.455923     28.47124   
    6       7  0.886736 -1.600869       0.986635    52.572016    15.594822   
    7       8  0.896672 -1.639343       0.992373     36.65884     9.454046   
    8       9  0.903531 -1.677747       0.994643    26.397061     6.024431   
    9      10  0.906045 -1.706474       0.995895    19.521431     3.975959   
    10     11  0.905582 -1.731076       0.994806    14.767291     2.623921   
    11     13  0.914482 -1.751347       0.997466     8.941254     1.205108   
    12     15  0.912684 -1.771227       0.994189     5.759987     0.568044   
    13     17  0.912188 -1.774908       0.990829     3.905403     0.273242   
    14     19  0.907649 -1.774186       0.989457     2.766824     0.135454   
    
             max(k)  
    0   5665.102661  
    1   3009.058821  
    2   1916.810605  
    3   1332.762771  
    4    984.036824  
    5    752.959999  
    6    591.514192  
    7    475.817182  
    8    389.237531  
    9    322.823838  
    10   270.867416  
    11   196.222414  
    12   146.575349  
    13   112.189052  
    14    87.594344  
    [92mSelected power to have scale free network is 9.[0m
    [96mcalculating adjacency matrix ...[0m
    	Done..
    
    [96mcalculating TOM similarity matrix ...[0m
    	Done..
    
    [96mGoing through the merge tree...[0m
    ..cutHeight not given, setting it to 0.996  ===>  99% of the (truncated) height range in dendro.
    	Done..
    
    [96mCalculating 34 module eigengenes in given set...[0m
    	Done..
    
    mergeCloseModules: Merging modules whose distance is less than 0.2
    fixDataStructure: data is not a Dictionary: converting it into one.
    multiSetMEs: Calculating module MEs.
      Working on set 1 ...
    [96mCalculating 34 module eigengenes in given set...[0m
    	Done..
    
    multiSetMEs: Calculating module MEs.
      Working on set 1 ...
    [96mCalculating 33 module eigengenes in given set...[0m
    	Done..
    
    multiSetMEs: Calculating module MEs.
      Working on set 1 ...
    [96mCalculating 31 module eigengenes in given set...[0m
    	Done..
    
    multiSetMEs: Calculating module MEs.
      Working on set 1 ...
    [96mCalculating 30 module eigengenes in given set...[0m
    	Done..
    
      Calculating new MEs...
    multiSetMEs: Calculating module MEs.
      Working on set 1 ...
    [96mCalculating 30 module eigengenes in given set...[0m
    	Done..
    
    [96mCalculating 30 module eigengenes in given set...[0m
    	Done..
    
    fixDataStructure: data is not a Dictionary: converting it into one.
    orderMEs: order not given, calculating using given set 0
    	Done running WGCNA..
    



    
![png](Quick_Start_files/Quick_Start_5_1.png)
    


We also can merge two previous steps by calling `runWGCNA()` function.

## Relating modules to external information and identifying important genes
PyWGCNA gather some important analysis after identifying modules in `analyseWGCNA()` function including:

1. Quantifying module–trait relationship 
2. Gene relationship to trait and modules

keep in your mind before start analysing don't forget to add any information you have about samples or genes.

For showing module relationship heatmap, PyWGCNA needs user to indicate color from [Matplotlib colors](https://matplotlib.org/stable/gallery/color/named_colors.html) for metadata by using `setMetadataColor()` function.

You also can select which data Trait in which order you wish to show in module eigen gene heatmap


```python
pyWGCNA_5xFAD.updateMetadata(path='5xFAD_paper/metaData', 
                            sep='\t')
# add color for metadata
pyWGCNA_5xFAD.setMetadataColor('Sex', {'Female': 'green',
                                       'Male': 'yellow'})
pyWGCNA_5xFAD.setMetadataColor('Genotype', {'5xFADWT': 'darkviolet',
                                            '5xFADHEMI': 'deeppink'})
pyWGCNA_5xFAD.setMetadataColor('Age', {'4mon': 'thistle',
                                       '8mon': 'plum',
                                       '12mon': 'violet',
                                       '18mon': 'purple'})
pyWGCNA_5xFAD.setMetadataColor('Tissue', {'Hippocampus': 'red',
                                          'Cortex': 'blue'})

geneList = PyWGCNA.getGeneList(dataset='mmusculus_gene_ensembl',
                               attributes=['ensembl_transcript_id', 
                                           'mgi_symbol', 
                                           'ensembl_gene_id', 
                                           'ensembl_peptide_id'])


pyWGCNA_5xFAD = pyWGCNA_5xFAD.analyseWGCNA(geneList=geneList)
```

    [1m[94mAnalysing WGCNA...[0m
    [96mCalculating module trait relationship ...[0m
    	Done..
    
    [96mAdding genes for each modules ...[0m
    	Done..
    
    [96mAdding gene name to gene module lists...[0m
    	Done..
    
    [96mplotting module heatmap eigengene...[0m
    	Done..
    
    [96mplotting Go term for each module...[0m


    2022-03-07 02:02:36,395 Warning: No enrich terms using library GO_Biological_Process_2021 when cutoff = 0.5
    2022-03-07 02:02:58,942 Warning: No enrich terms using library GO_Biological_Process_2021 when cutoff = 0.5


    	Done..
    


## Saving and loading your PyWGCNA
you can save or load your PyWGCNA object with `saveWGCNA()` or `readWGCNA()` function.


```python
pyWGCNA_5xFAD.saveWGCNA()
```

    [1m[94mSaving WGCNA as 5xFAD.p[0m


you can also load your PyWGCNA object with `readWGCNA()` function.


```python
import PyWGCNA
pyWGCNA_5xFAD = PyWGCNA.readWGCNA("5xFAD.p")
```

    [1m[94mReading WGCNA done![0m

