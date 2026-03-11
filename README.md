1️⃣ GitHub README
Project Title

Single-Cell RNA-seq Analysis of Human Pancreatic Cells Using Seurat

Overview

This project performs single-cell RNA sequencing (scRNA-seq) analysis of human pancreatic cells to identify transcriptionally distinct cell populations and explore gene expression variability across cell types.

The analysis includes:

Quality control filtering

Gene expression normalization

Identification of highly variable genes

Dimensionality reduction (PCA)

Graph-based clustering

UMAP visualization of cell populations

Marker gene discovery

The workflow was implemented using the Seurat package in R.

Dataset

The dataset used is the panc8 pancreatic single-cell dataset, accessed through the SeuratData library.

Dataset summary:

Feature	Value
Cells analyzed	~10,800
Genes	~34,000
Cell types	8
Species	Human

Major cell populations include:

acinar

alpha

beta

delta

ductal

endothelial

activated stellate

gamma

Workflow
1 Data Loading

The pancreatic dataset was loaded from the SeuratData repository and converted to the latest Seurat object format.

2 Quality Control

Cells were filtered based on:

nFeature_RNA > 200
nFeature_RNA < 5000
percent.mt < 5

This removes low-quality or damaged cells.

QC Visualization

These plots show:

gene counts per cell

RNA molecule counts

mitochondrial gene percentage

Variable Gene Identification

Highly variable genes were identified using variance stabilization.

Red points represent genes with high expression variance across cells.

Example genes identified:

SPP1

OLFM4

CFTR

MMP7

RBP4

Dimensionality Reduction

Principal Component Analysis (PCA) was used to reduce dimensionality.

Elbow Plot

This plot determines the number of principal components to retain.

PCA Heatmap

The heatmap shows genes contributing most strongly to each principal component.

Clustering

Cells were grouped using a graph-based clustering algorithm.

Result:

18 clusters identified

UMAP Visualization

UMAP projects high-dimensional gene expression into a 2-dimensional space.

Cluster Visualization

Each color represents a transcriptionally distinct cluster.

Known Cell Types

Clusters correspond closely to known pancreatic cell populations.

Key Findings

1️⃣ Distinct pancreatic cell populations were clearly separated in UMAP space.

2️⃣ Highly variable genes strongly contributed to cluster separation.

3️⃣ Known biological cell types aligned well with computational clusters.

Technologies Used
Tool	Purpose
R	Statistical computing
Seurat	Single-cell RNA analysis
SeuratData	Dataset access
dplyr	Data manipulation
ggplot2	Visualization
Project Structure
project/
│
├── scripts/
│   analysis.R
│
├── figures/
│   qc_violin.png
│   variable_genes.png
│   elbow_plot.png
│   pca_heatmap.png
│   umap_clusters.png
│   umap_celltypes.png
│
└── pancreas_cluster_markers.csv
