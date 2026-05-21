# Doublet Detection Benchmark for hm-6k scRNA-seq Data

This project compares three doublet detection methods on the human-mouse mixture single-cell RNA-seq dataset:

- scDblFinder
- DoubletFinder
- Scrublet

The aim of this project is to evaluate how different doublet detection tools perform on the same annotated dataset and to examine the effect of doublet removal on downstream single-cell analysis.

R v4.5.2 was used for Seurat preprocessing, scDblFinder, DoubletFinder, benchmarking, and plotting.  
Python v3.9.25 was used for running Scrublet.

## Dataset
This project uses the hm-6k human-mouse mixture scRNA-seq dataset.

## Workflow

The analysis includes:

1. Loading the hm-6k dataset

2. Creating a Seurat object

3. Quality control and filtering

4. Normalisation， feature selection，PCA, clustering, and UMAP

5. Running scDblFinder

6. Running DoubletFinder

7. Importing Scrublet results

8. Comparing doublet detection performance

9. Removing predicted doublets

10. Re-running downstream clustering analysis

11.Saving figures and result tables
