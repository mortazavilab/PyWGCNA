# Data input, cleaning and pre-processing

This is the first step of any network analysis. 
We show here how to load typical expression data, 
pre-process them into a format suitable for 
network analysis, and clean the data by removing 
obvious outlier samples as well as genes and 
samples with excessive numbers of missing entries.

* [Data Input](Data_format.md#data-input)
  * [anndata format](Data_format.md#expression-data-gene-and-sample-information-all-together-in-anndata-format)
  * [separate format](Data_format.md#expression-data-gene-and-sample-information-separately)
    * [Gene Expression](Data_format.md#gene-expression)
    * [Gene Information](Data_format.md#gene-information)
    * [Sample Information](Data_format.md#sample-information)
    * [Other parameters](Data_format.md#other-parameters)
* [Data cleaning and pre-processing](Data_format.md#data-cleaning-and-pre-processing)

## Data Input
We store **raw** expression data along information in [anndata](https://anndata.readthedocs.io/en/latest/) format in `geneExpr` variable.
you can pass your expression data, gene and sample information all together or separately:

### expression data, gene and sample information all together in anndata format
If you already have your expression data in anndata format you can define your pyWGCNA object by passing your variable in `anndata` format. 
keep in mind X should be expression matrix. var is gene information and obs is sample information.

### expression data, gene and sample information separately
you can pass the paths that store each information or the table contains them.

#### Gene Expression
The expression data is a table which the rows are samples and columns are genes.
The first column (index of dataframe) is going to be sample id or sample name and first column (column of dataframe) should be gene id or gene name which both of them should be unique.

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>sample_id</th>
      <th>ENSMUSG00000000003</th>
      <th>ENSMUSG00000000028</th>
      <th>ENSMUSG00000000031</th>
      <th>ENSMUSG00000000037</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>sample_11615</td>
      <th>12.04</th>
      <th>11.56</th>
      <th>16.06</th>
      <th>13.18</th>
    </tr>
    <tr>
      <td>sample_11616</td>
      <th>1.35</th>
      <th>1.63</th>
      <th>1.28</th>
      <th>1</th>
    </tr>
  </tbody>
</table>
</div>

#### Gene Information
The gene information is a table which contains additional information about each genes. 
First column should be your index which should be the same name as first column of gene expression data (gene ID).


<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>gene_id</th>
      <th>gene_name</th>
      <th>gene_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>ENSMUSG00000000003</td>
      <th>Pbsn</th>
      <th>protein_coding</th>
    </tr>
    <tr>
      <td>ENSMUSG00000000028</td>
      <th>Cdc45</th>
      <th>protein_coding</th>
    </tr>
    <tr>
      <td>ENSMUSG00000000031</td>
      <th>H19</th>
      <th>lncRNA</th>
    </tr>
    <tr>
      <td>ENSMUSG00000000037</td>
      <th>Scml2</th>
      <th>protein_coding</th>
    </tr>
  </tbody>
</table>
</div>

#### Sample Information
The sample information is a table which contains additional information about each sample. 
First column should be your index which should be the same name as first row of gene expression data (sample ID).

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Sample_id</th>
      <th>Age</th>
      <th>Tissue</th>
      <th>Sex</th>
      <th>Genotype</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>sample_11615</td>
      <td>4mon</td>
      <td>Cortex</td>
      <td>Female</td>
      <td>5xFADHEMI</td>
    </tr>
    <tr>
      <td>sample_11616</td>
      <td>4mon</td>
      <td>Cortex</td>
      <td>Female</td>
      <td>5xFADWT</td>
    </tr>
  </tbody>
</table>
</div>

### Other parameters
These are other parameters we suggest checking them before starting any analysis. 
* **name**: name of the WGCNA we used to visualize data (default: `WGCNA`)

* **save**: define whether you want to save result of important steps or not (If you want to set it 
`TRUE` you should have a write access on the output directory)

* **outputPath**: define where you want to save your data, otherwise it will be store near the code. 

* **TPMcutoff**: cut off for removing genes that expressed under this number along samples

* **networkType** : Type of networks (default: `signed hybrid` and Options: `unsigned`, `signed` and `signed hybrid`)

* **adjacencyType**: Type of adjacency matrix (default: `signed hybrid` and Options: `unsigned`, `signed` and `signed hybrid`)


* **TOMType**: Type of topological overlap matrix(TOM) (default: `signed` and Options: `unsigned` and `signed`)

For depth-in documents look at [here](https://mortazavilab.github.io/PyWGCNA/html/PyWGCNA.html).

## Data cleaning and pre-processing

PyWGCNA checks data for genes and samples with too many missing values.
1. Remove genes without any expression more than `TPMcutoff` value (default one) across all samples.
2. `goodSamplesGenes()` function to find genes and samples with too many missing values.
3. Cluster the samples (use [Hierarchical clustering](https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html#module-scipy.cluster.hierarchy)
from [scipy](https://scipy.org/)) to see if there are any obvious outliers. 
you can define value the height by `cut` value. By default, we don't remove any sample by hierarchical clustering
