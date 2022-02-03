# PyWGCNA

PyWGCNA is a Python library designed to do Weighted correlation network analysis (WGCNA) 
can be used for finding clusters (modules) of highly correlated genes, for summarizing 
such clusters using the module eigengene for relating modules to one another and 
to external sample traits (using eigengene network methodology), and for calculating 
module membership measures. Users can also compare WGCNA from different datasets
including single cell gene markers.

## Installation

To install PyWGCNA, python version 3.7 or greater is required.

### Install from PyPi
Install the most recent release, run

`pip install PyWGCNA`

### Install with the most recent commits
git cloning the [PyWGCNA repository](https://github.com/mortazavilab/PyWGCNA), going to the PyWGCNA directory, run

`pip install .`

## Tutorials

- [Data input, cleaning and pre-processing](tutorials/Data%20format.md)
- [Quick Start](tutorials/Quick%20Start.ipynb): How to load data into PyWGCNA and find modules and analyse them

