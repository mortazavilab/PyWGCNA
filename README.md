# PyWGCNA

PyWGCNA is a Python library designed to do weighted correlation network analysis (WGCNA). It
can be used for finding clusters (modules) of highly correlated genes, for summarizing
such clusters using the module eigengene, for relating modules to one another and
to external sample traits (using eigengene network methodology), and for calculating
module membership measures. Users can also compare WGCNA networks from different datasets, or to
external gene lists, to assess the conservation or functional enrichment of each module.

![PyWGCNA overview](docs/PyWGCNA_overview.png)

## Documentation
PyWGCNA's full documentation can be found [here](https://mortazavilab.github.io/PyWGCNA/)

## Installation

To install PyWGCNA, Python version 3.7 or greater is required.

### Install from PyPi (recommended)
To install the most recent release, run

`pip install PyWGCNA`

### Install with the most recent commits
* Git clone the [PyWGCNA repository](https://github.com/mortazavilab/PyWGCNA), cd to the `PyWGCNA` directory, and run

`pip install .`

## Tutorials

- [Data input, cleaning and pre-processing](tutorials/Data_format.md): How to format, clean and preprocess your input data for PyWGCNA
- [Quick Start](tutorials/Quick_Start.ipynb): How to load data into PyWGCNA, find modules, and analyze them
- [PyWGCNA object](tutorials/PyWGCNA_object.ipynb): How to interact with PyWGCNA objects and some parameters we have them in the object and how you can access them
- [Compare two PyWGCNA objects](tutorials/Comparison_two_PyWGCNAs.ipynb): How to compare two PyWGCNA objects
- [Compare more than two PyWGCNA objects](tutorials/Comparison_multi_PyWGCNAs.ipynb): How to compare three PyWGCNA objects
- [Compare PyWGCNA objects to gene marker list](tutorials/Comparison_PyWGCNA_geneMarker.ipynb): How to compare PyWGCNA objects to external gene lists (here shown on marker genes from single-cell data)
- [Functional enrichment analysis ](tutorials/functional_enrichment_analysis.ipynb): How to perform functional enrichment analysis using databases such as GO, KEGG, and REACTOME in PyWGCNA object
- [Visualize modules as network](tutorials/network_analysis.ipynb): How to visualize PyWGCNA objects as a network
- [Recover Protein-Protein Interaction](tutorials/protein_protein_interaction.ipynb): How to find and plot PPI using STRING database.

## Suggested Reading

If you are unfamiliar with R refrence WGCNA, we suggest reading the original WGCNA publication:

- [WGCNA: an R package for weighted correlation network analysis](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)
- [Describe various features of the WGCNA package for the R programming language](https://peterlangfelder.com/)

## Cite

PyWGCNA is now online in Bioinformatics. Please cite our paper when using PyWGCNA:

Narges Rezaie and others, PyWGCNA: A Python package for weighted gene co-expression network analysis, Bioinformatics, 2023;
[https://doi.org/10.1093/bioinformatics/btad415](https://doi.org/10.1093/bioinformatics/btad415)
