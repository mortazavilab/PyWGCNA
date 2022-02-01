# Data input, cleaning and pre-processing

This is the first step of any network analysis. 
We show here how to load typical expression data, 
pre-process them into a format suitable for 
network analysis, and clean the data by removing 
obvious outlier samples as well as genes and 
samples with excessive numbers of missing entries.

* [Data Input](Data%20format.md#data-input)
    - [Gene Expression](Data%20format.md#gene-expression)
    - [Gene Information](Data%20format.md#gene-information)
    - [Sample Information](Data%20format.md#sample-information)
    - [Other parameters](Data%20format.md#other-parameters)
* [Data cleaning and pre-processing](Data%20format.md#data-input-cleaning-and-pre-processing)

## Data Input

Here we explain important data we used as an 
input.

###Gene Expression
The expression data is a table which the rows are 
genes and columns are samples, the first column is 
gonna be gene_id or gene name and
first column should be sample id.

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>gene_id</th>
      <th>sample_11615</th>
      <th>sample_11616</th>
      <th>sample_11617</th>
      <th>sample_11618</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>ENSMUSG00000000001.4</td>
      <th>12.04</th>
      <th>11.56</th>
      <th>16.06</th>
      <th>13.18</th>
    </tr>
    <tr>
      <td>ENSMUSG00000000028.15</td>
      <th>1.35</th>
      <th>1.63</th>
      <th>1.28</th>
      <th>1</th>
    </tr>
  </tbody>
</table>
</div>

### Gene Information
The gene information is a table which contains 
additional information about each genes. It should 
have a same order as gene expression matrix. 


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
      <td>ENSMUSG00000000001.4</td>
      <th>Gnai3</th>
      <th>Protein coding</th>
    </tr>
    <tr>
      <td>ENSMUSG00000000028.15</td>
      <th>Cdc45</th>
      <th>Protein coding</th>
    </tr>
  </tbody>
</table>
</div>

### Sample Information
The sample information is a table which contains 
additional information about each samples. It should 
have a same order as gene expression matrix.

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Sample_id</th>
      <th>Sex</th>
      <th>Age</th>
      <th>Genotype</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>sample_11615</td>
      <td>Female</td>
      <td>4 month</td>
      <td>5xFAD</td>
    </tr>
    <tr>
      <td>sample_11616</td>
      <td>Male</td>
      <td>4 month</td>
      <td>5xFAD</td>
    </tr>
    <tr>
      <td>sample_11617</td>
      <td>Female</td>
      <td>4 month</td>
      <td>BL6</td>
    </tr>
    <tr>
      <td>sample_11618</td>
      <td>Male</td>
      <td>4 month</td>
      <td>BL6</td>
    </tr>
  </tbody>
</table>
</div>

### Other parameters
* **name**: name of the WGCNA we used to visualize
data (default: 'WGCNA')
* **save**: define whether you want to save result 
of important steps or not (If you want to set it 
TRUE you should have a write access on the output 
directory)
* **outputPath**: define where you want to save 
your data, otherwise it will be store near the 
code. 
* **networkType** : Type of networks (Options: 
"unsigned", "signed" and "signed hybrid")

* **adjacencyType**: Type of adjacency matrix 
(Options: "unsigned", "signed" and "signed hybrid")

* **TOMType**: Type of topological overlap matrix
(TOM) (Options: "NA", "unsigned", "signed")



## Data cleaning and pre-processing

PyWGCNA checks data for genes and samples 
with too many missing values.
1. Remove genes without any expression more 
than one (or you can define the number by 
changing `TPMcutoff` value) across all samples.
2. `goodSamplesGenes()` function to find 
genes and samples with too many missing values.
3. Cluster the samples (use [Hierarchical clustering](https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html#module-scipy.cluster.hierarchy)
from [scipy](https://scipy.org/)) to see if 
there are any obvious outliers.