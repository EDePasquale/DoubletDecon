# DoubletDecon #

A cell-state aware tool for removing doublets from single-cell RNA-seq data

![logo](http://www.altanalyze.org/DoubletDecon/wordcloud.png)

See our [bioRxiv](https://www.biorxiv.org/content/early/2018/07/08/364810) for more information on DoubletDecon.


# Installation

```javascript
devtools::install_github('EDePasquale/DoubletDecon')
```

# Dependencies # 

DoubletDecon requires the following R packages:

* DeconRNASeq
* gplots
* dplyr
* MCL
* clusterProfiler
* mygene
* hopach
* as.color

Additionally, the use of the cell cycle removal option requires an internet connection.


# Usage #

## Seurat data only: ##

```javascript
Seurat_Pre_Process(expressionFile, genesFile, clustersFile)
```

#### Arguments ####

* expressionFile: Normalized expression matrix as a .txt file from Seurat's NormalizeData() function
* genesFile: Top marker gene list as a .txt file from Seurat's top_n() function
* clustersFile: Cluster identities as a .txt file from Seurat object @ident 


#### Value ####

* newExpressionFile - Seurat expression file in ICGS format (used as 'rawDataFile')
* newGroupsFile - Groups file ICGS format (used as 'groupsFile')

## Seurat and ICGS data: ##

```javascript
Main_Doublet_Decon(rawDataFile, groupsFile, filename, location,
  fullDataFile = NULL, removeCC = FALSE, species = "mmu", rhop = 1,
  write = TRUE, recluster = "doublets_decon", PMF = TRUE,
  useFull = FALSE, heatmap = TRUE)
```

#### Arguments ####

* rawDataFile: Name of file containing ICGS or Seurat expression data (gene by cell)
* groupsFile: Name of file containing group assignments (3 column: cell, group(numeric), group(numeric or character))
* filename: Unique filename to be incorporated into the names of outputs from the functions.
* location: Directory where output should be stored
* fullDataFile: Name of file containing full expression data (gene by cell). Default is NULL.
* removeCC: Remove cell cycle gene cluster by KEGG enrichment. Default is FALSE.
* species: Species as scientific species name, KEGG ID, three letter species abbreviation, or NCBI ID. Default is "mmu".
* rhop: x in mean+x*SD to determine upper cutoff for correlation in the blacklist. Default is 1.
* write: Write output files as .txt files. Default is TRUE.
* recluster: Recluster deconvolution classified doublets and non-doublets seperately using hopach or deconvolution classifications.
* PMF: Use step 3 (unique gene expression) in doublet determination criteria. Default is TRUE.
* useFull: Use full gene list for PMF analysis. Requires fullDataFile. Default is FALSE.
* heatmap: Boolean value for whether to generate heatmaps. Default is TRUE. Can be slow to datasets larger than ~3000 cells.

#### Value ####

* data_processed = new expression file (cleaned).
* groups_processed = new groups file (cleaned).
* PMF_results = pseudo marker finder t-test results (gene by cluster).
* DRS_doublet_table = each cell and whether it is called a doublet by deconvolution analysis.
* DRS_results = results of deconvolution analysis (cell by cluster) in percentages.
* Decon_called_freq = percentage of doublets called in each cluster by deconvolution analysis.
* Final_doublets_groups = new groups file containing only doublets.
* Final_nondoublets_groups = new groups file containing only non doublets.
