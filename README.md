# DoubletDecon #

A cell-state aware tool for removing doublets from single-cell RNA-seq data

![logo](http://www.altanalyze.org/DoubletDecon/wordcloud.png)

See our [bioRxiv](https://www.biorxiv.org/content/early/2018/07/08/364810) for more information on DoubletDecon.

# Updates - Version 1.0.2 : January 9th, 2019 #
 * General bug fixes affecting final groups file and final expression file output.
 * Added user option to specify minimum number of unique genes to Rescue a putative doublet cluster (previously set at 4).
 * Speed up run time for users who do not use the Rescue step (PMF=FALSE).
 * Remove requirement for 'as.color' function.

## Version 1.0.1 : December 26th, 2018
 * Additional "Remove" step option to create synthetic doublet centroids with 30%/70% and 70%/30% parent cell contribution instead of simply 50%/50% (only50=FALSE).
 * Heatmap generation corrected for large datasets (>5000 cells).
 * "Rescue" step modification from t-tests for all clusters to ANOVA with Tukey post-hoc test in only putative doublet clusters. Minimum of 4 unique genes as hardcoded default.
 * "Rescue" step now allows for sampling of clusters evenly or proportional to cluster size when using the full expression matrix.
 * Hopach removed as a "Recluster" option; does not work with improved "Rescue" step. Subsequently removed the DeconCalledFreq table as written and returned output.
 * Log file is generated with a unique ID for each run of DoubletDecon.
 * Catch error in mcl() function and quits DoubletDecon with warning to choose a different rhop value.
 * Synthetic doublet deconvolution values output for quality control (Synth_doublet_info)


# Installation #

Run the following code to install the package using devtools:

```
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
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
 
These can be installed with:

```
source("https://bioconductor.org/biocLite.R")
biocLite(c("DeconRNASeq", "clusterProfiler", "hopach", "mygene"))
install.packages("MCL")
```

Additionally, the use of the cell cycle removal option requires an internet connection.
 
 
# Usage #

## Seurat data only: ##

```javascript
Seurat_Pre_Process(expressionFile, genesFile, clustersFile)
```

#### Arguments ####

* expressionFile: Normalized expression matrix or counts file as a .txt file (expression from Seurat's NormalizeData() function)
* genesFile: Top marker gene list as a .txt file from Seurat's top_n() function
* clustersFile: Cluster identities as a .txt file from Seurat object @ident

#### Value ####

* newExpressionFile - Seurat expression file in ICGS format (used as 'rawDataFile')
* newGroupsFile - Groups file ICGS format (used as 'groupsFile')

## Seurat and ICGS data: ##

```javascript
Main_Doublet_Decon(rawDataFile, groupsFile, filename, location,
  fullDataFile = NULL, removeCC = FALSE, species = "mmu", rhop = 1,
  write = TRUE, PMF = TRUE, useFull = FALSE, heatmap = TRUE, centroids=FALSE, num_doubs=30, 
  downsample="none", sample_num=NULL, only50=TRUE, min_uniq=4)
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
* centroids: Use centroids as references in deconvolution instead of the default medoids.
* num_doubs: The user defined number of doublets to make for each pair of clusters. Default is 30.
* downsample: allows for downsampling of cells when using full expression matrix (use with large datasets), default is "none".
* sample_num: number of cells per cluster with downsampling with "even", percent of cluster with "prop".
* only50: use only synthetic doublets created with 50%/50% mix of parent cells, as opposed to the extended option of 30%/70% and 70%/30%, default is TRUE.
* min_uniq: minimum number of unique genes required for a cluster to be rescued, default is 4.

#### Value ####

* data_processed = new expression file (cleaned).
* groups_processed = new groups file (cleaned).
* PMF_results = pseudo marker finder t-test results (gene by cluster).
* DRS_doublet_table = each cell and whether it is called a doublet by deconvolution analysis.
* DRS_results = results of deconvolution analysis (cell by cluster) in percentages.
* Decon_called_freq = percentage of doublets called in each cluster by deconvolution analysis.
* Final_doublets_groups = new groups file containing only doublets.
* Final_nondoublets_groups = new groups file containing only non doublets.
* Synth_doublet_info = synthetic doublet deconvolution values output for quality control.


# Example #

Data for this example can be found in this GitHub, and in combination with the below function calls, can reproduce the results from Figure 5 (Identification of Experimentally Verified Doublets from PBMC) of the bioRxiv pre-print.

Update: This figure may change for the final publication as some of the functions work differently in Version 1.0.1. The example will still work but it will no longer exactly match the preprint.

```javascript
location="/Users/xxx/xxx/" #Update as needed 
expressionFile=paste0(location, "counts.txt")
genesFile=paste0(location, "Top50Genes.txt")
clustersFile=paste0(location, "Cluster.txt")

newFiles=Seurat_Pre_Process(expressionFile, genesFile, clustersFile)

results=Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                           groupsFile=newFiles$newGroupsFile, 
                           filename="PBMC_example", 
                           location=location,
                           fullDataFile=NULL, 
                           removeCC=FALSE, 
                           species="hsa", 
                           rhop=1.1, 
                           write=TRUE, 
                           PMF=TRUE, 
                           useFull=FALSE, 
                           heatmap=FALSE,
                           centroids=TRUE,
                           num_doubs=100, 
                           downsample="none",
                           sample_num=NULL,
                           only50=TRUE,
                           min_uniq=4)
```                           
