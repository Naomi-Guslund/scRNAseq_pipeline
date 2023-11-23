# scRNAseq_pipeline
Taking raw sequencing files and taking them through a pipeline until visualising cell clusters

# Characterising the splenic cell populations in the Atlantic cod using scRNAseq
Here I provide example scripts used to analyse scRNAseq data in order to characterise splenic cells of the Atlantic cod, a non-model organism.
Order of scripts provided:
1) From_fastq_to_Digital_expression_matrix - bash script
2) LoadingData.html - to view, download the file and open in your preferred directory
3) PreprocessingWorkflow_DimensionReduction.html -to view, download the file and open in your preferred directory

## Biological background
Atlantic Cod (Gadus morhua) has lost the major histocompatibility complex class II presentation pathway. Identification and characterization of immune cell subsets is needed to understand how this alternative immune system functions. Here, we use single-cell RNA sequencing to examine the cellular heterogeneity in Atlantic cod spleen. We analysed spleens from 34 Atlantic cod at 12 timepoints during a vaccination and immune challenge study of V. anguillarum, thus capturing a broad immune status. The 3 files demonstrate how raw sequencing files are demultiplexed and aligned to the genome (From_fastq_to_Digital_expression_matrix), how data is loaded into R and organised into a seurat object (LoadingData.html). The final R script demonstrates how normalisation is carried out, how filtering of low quality cells are removed, and how cells can be visualised on a UMAP (PreprocessingWorkflow_DimensionReduction.html).

## Dataset
The sequencing data is available at the ENA repository with Accession number PRJEB47815.
The project analyses a dataset of spleen cells from the Atlantic cod, captured using the Dolomite Single-Cell RNA-Seq System, sequenced on the Illumina NextSeq 500. The Drop-seq sequencing libraries produce paired-end reads: read 1 contains both a cell barcode and a molecular barcode (also known as a UMI); read 2 is aligned to the reference genome- the most recent version of the Atlantic cod genome, gadMor3 (RefSeq accession GCF_902167405.1) using STAR alignment. The raw sequence data was demultiplexed by folowing the "Drop-seq Core Computational Protocol" written by James Nemesh (version 2.0.0 (9/28/18)), producing a “digital expression matrix” containing integer counts of the number of transcripts for each gene (by row), in each cell (cell barcodes as columns).

### Analysis steps in R
The resulting count matrices for each sample were loaded into Seurat (version 4.0.2). The count matrices from each of the 34 cod were integrated into one Seurat object through anchor identification using the ‘FindIntergationAnchors’ and ‘IntergrateData’ Seurat functions. After filtering away low-quality cells and cell multiplets, we derived a gene expression matrix of 19,279 genes across 56,994 splenic cells. Data were normalized using the ‘NormalizeData’ function, variable features were identified using the ‘FindVariableFeatures’ function, and features were scaled using the ‘ScaleData’ function. PCA analysis was performed, and the most significant principal components were identified.
 
Unsupervised dimension reduction was performed using Uniform Manifold Approximation and Projection (UMAP) the cell clusters were identified using the ‘FindNeighbors’ and ‘FindClusters’ functions. Cluster identities can be assigned by assessing differential expression of marker genes between different clusters using the biomarkers detected from the ‘FindAllMarkers’ function.

## Published paper
Guslund, N.C., Krabberød, A.K., Nørstebø, S.F. et al. Lymphocyte subsets in Atlantic cod (Gadus morhua) interrogated by single-cell sequencing. Commun Biol 5, 689 (2022). https://doi.org/10.1038/s42003-022-03645-w

## Resources that were used to analyse data 
Drop-seq Core Computational Protocol - https://raw.githubusercontent.com/broadinstitute/Drop-seq/v2.4.0/doc/Drop-seq_Alignment_Cookbook.pdf
Online tutorial for handling scRNAseq data - https://satijalab.org/seurat/
Online tutorial for trajectory inference workflow - https://dynverse.org/users/2-quick_start/
