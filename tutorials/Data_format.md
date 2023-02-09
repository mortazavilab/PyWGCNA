# Data input, cleaning and pre-processing

This is the first step of any network analysis.
We show here how to load typical expression data,
pre-process them into a format suitable for
network analysis, and clean the data by removing
obvious outlier samples and genes genes.

* [Data Input](#input_data)
  * [AnnData format](#anndata)
  * [Separate matrices](#separate_mat)
    * [Gene expression](#gene_exp)
    * [Gene metadata](#gene_meta)
    * [Sample metadata](#sample_meta)
    * [Other parameters](#params)
* [Data cleaning and pre-processing](#data_preproc)

## <a name="input_data"></a>Input data format

We store **raw** expression data along information in [AnnData](https://anndata.readthedocs.io/en/latest/) format in the `geneExpr` variable. Gene expression data, gene metadata, and sample metadata can either be passed to PyWGCNA all together in an AnnData object, or separately as a series of matrices.

### <a name="anndata"></a>AnnData format
If you already have your expression data in AnnData format you can define your PyWGCNA object by passing your variable in `AnnData` format.
Keep in mind `AnnData.X` should be the expression matrix, `AnnData.var` should contain information about each gene, and `AnnData.obs` should contain information about each sample. You can read more about the AnnData format [here](https://anndata.readthedocs.io/en/latest/)

### <a name="separate_mat"></a>Separate matrices for gene expression, sample metadata, and gene metadata

The user can pass individual file paths for gene expression, sample metadata, and gene metadata, in the formats specified below.

#### <a name="gene_exp"></a>Gene expression
The expression table should be formatted such that the rows correspond to samples and the columns correspond to genes.
The first column should represent the sample id or sample name. The following columns should contain gene ids or gene names which are all unique.

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

#### <a name="gene_meta"></a>Gene metadata

The gene metadata is a table which contains additional information about each gene, such as gene biotype or gene length.
Each row should represent a gene and each column should represent a gene feature, where the first columns contains the same gene identifier that was used in the gene expression matrix
The rows should be in the same order as the columns of the gene expression matrix, or
the user can specify `order=False`.

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

#### <a name="sample_meta"></a>Sample metadata

The sample metadata is a table which contains additional information about each sample, such as timepoint or genotype.
Each row should represent a sample and each column should represent a metadata feature, where the first columns contains the same sample identifier that was used in the gene expression matrix
The rows should be in the same order as the rows of the gene expression matrix, or
the user can specify `order=False`.

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

### <a name="params"></a>Other parameters
These are other parameters that can be specified.

* **name**: Name of the WGCNA used to visualize data (default: `WGCNA`)

* **save**: Whether to save the results of important steps or not (If you want to set it
`True` you should have a write access on the output directory)

* **outputPath**: Where to save your data, otherwise it will be stored in the same directory as the code.

* **TPMcutoff**: TPM cutoff for removing genes

* **networkType** : Type of network to generate ({`unsigned`, `signed` and `signed hybrid`}, default: `signed hybrid`)

* **adjacencyType**: Type of adjacency matrix to use ({`unsigned`, `signed` and `signed hybrid`}, default: `signed hybrid`)

* **TOMType**: Type of topological overlap matrix(TOM) to use ({`unsigned`, `signed`}, default: `signed`)

For depth-in documentation on these parameters see [here](https://mortazavilab.github.io/PyWGCNA/html/PyWGCNA.html).

## <a name="data_preproc"></a>Data cleaning and preprocessing

PyWGCNA can clean the input data according to the following criteria:
1. Remove genes without any expression more than `TPMcutoff` value (default one) across all samples.
2. Find genes and samples `goodSamplesGenes()` function to find genes and samples with too many missing values.
3. Cluster the samples (uses [hierarchical clustering](https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html#module-scipy.cluster.hierarchy)
from [scipy](https://scipy.org/)) to see if there are any obvious outliers. The user can define value the height by specifying the `cut` value. By default, no samples are removed by hierarchical clustering
