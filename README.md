# DoubletDecon #

Deconvoluting doublets from single-cell RNA-sequencing data

![logo](http://www.altanalyze.org/DoubletDecon/wordcloud.png)

See our [Cell Reports paper](https://www.cell.com/cell-reports/fulltext/S2211-1247(19)31286-0) for more information on DoubletDecon. Also see our [bioRxiv](https://www.biorxiv.org/content/early/2018/07/08/364810) for an older description of the algorithm.

NEW! See our protocol on [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.04.23.058156v1) for more description on how to use DoubletDecon.

# Updates: November 30th, 2020 #
 * Improved_Seurat_Pre_Process() — This function now automatically upgrades your Seurat object to work with the version of Seurat you have installed (this is to help with the transition to Seurat version 4). Users are now able to specify if the expression values should be pulled from the @counts, @data, or @scaled.data slot in the Seurat object. Redundant code was removed.  The section of code using FindAllMarkers() has been updated to look for the column titled "avg_log2FC" in Seurat 4 and "avg_logFC" in Seurat 3 (3.2.2 and 3.9.9 tested). Cluster numbers should now not increment by 2, but by 1 as was originally intended.

# URGENT NOTE : July 2nd, 2020 #
  * It has been brought to my attention that the UI version of DoubletDecon distributed with the latest version (1.1.5) does not work. We are actively working to solve this as soon as possible. However, the command line version of DoubletDecon is still working. If you would like to use the UI version still, you can open R and use the command:
  
```
shiny::runGist('a81cdc2aea5742c08e5fc3fa66d47698', launch.browser=TRUE)
```
   This temporary solution will work until the application is fixed. Thank you for your patience!

# Updates - Version 1.1.5 : May 27th, 2020 #
  * Signed and notarized DoubletDeconUI app to help limit Gatekeeper issues with new MacOS
  * NEW! Vignette in the Wiki portion of this repository for playing with DoubletDecon on real data
  * Updated DESCRIPTION and Main_Doublet_Decon() function to improve the installation experience
  
# Updates - Version 1.1.4 : January 6th, 2020 #
  * Fixed bug in Improved_Seurat_Pre_Process caused by an incorrect assumption that cell names were in the first column and not the column names in the Seurat expression object
  * Added new parameter to Main_Doublet_Decon to allow for manual override of the automatic cores detection used in the 'rescue' step. The default is set to -1, which triggers automatic detection and should replicate the existing experience.

# Updates - Version 1.1.3 : November 6th, 2019 #
  * NEW! Integrated ICGS2_to_ICGS1() is now available to support input files from ICGS version 2. You should not have to make any changes in your DoubletDecon workflow to use ICGS version 2 instead of ICGS version 1
  * Improved_Seurat_Pre_Process() now loads dplyr at the beginning of the fuction (thanks chansigit for the feedback!)
  * Fixed bugs in Remove_Cell_Cycle() that was preventing it from running on certain datasets
  
# Updates - Version 1.1.2 : September 5th, 2019 #
  * NEW! Improved_Seurat_Pre_Process() is now available to replace Seurat_Pre_Process() for those who would prefer to work directly with a Seurat Object as input instead of individual files saved from a Seurat workflow. Workflows following the protocol found at https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html, from the provided script (seurat-3.0.R), or similar will be sufficent for this new function.
  * Resolved compatibility issues with Seurat version 3

# Updates - Version 1.1.1 : May 29th, 2019 #
 * Change default for only50 to FALSE (from TRUE) to reflect best practices in running DoubletDecon

# Updates - Version 1.1.0 : March 26th, 2019 #
 * NEW! DoubletDecon UI available in the GitHub repository (requires R3.5.0 or later and RStudio with 'shiny' package installed)
 * NEW! Improved Rescue step, Pseudo_Marker_Finder, now uses parallel processing and data chunking to improve speed and memory efficency. Results remain the same with the exception of no longer saving p-values (future release)
 * Remove downsample and sample_num
 * Upped num_doubs default value to 100 (from 30)
 * Require 'tidyr', 'R.utils', 'forrach', 'doParallel', 'stringr'. No longer require 'hopach'
 * Changed log_name_file to the value for filename, for compatibility with Windows operating systems
 * Automatically write processed data and groups files from Clean_Up_Data (even if write=FALSE) for use with new Rescue step (Pseudo_Marker_Finder)
 * Finalize switch to more granular doublet calls in the case of no Rescue step
 * Change name of cluster merging plot to "Cluster Merge" (from "Blacklist")
 * Fixed bug in Remove step when number of clusters equals 2 (Euclidean distance is used in the place of Pearson correlation)

# Updates - Version 1.0.2 : January 9th, 2019 #
 * General bug fixes affecting final groups file and final expression file output.
 * Added user option to specify minimum number of unique genes to Rescue a putative doublet cluster (previously set at 4).
 * Speed up run time for users who do not use the Rescue step (PMF=FALSE).
 * Remove requirement for 'as.color' function.

# Updates - Version 1.0.1 : December 26th, 2018 #
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
 * tidyr
 * R.utils
 * foreach
 * doParallel
 * stringr
 * Seurat (for Improved_Seurat_Pre_Process only)
 
These can be installed with:

```
source("https://bioconductor.org/biocLite.R")
biocLite(c("DeconRNASeq", "clusterProfiler", "hopach", "mygene", "tidyr", "R.utils", "foreach", "doParallel", "stringr"))
install.packages("MCL")
```

Additionally, the use of the cell cycle removal option requires an internet connection.
 
 
# Usage #

## Seurat data only: ##

```javascript
Improved_Seurat_Pre_Process(seuratObject, num_genes=50, write_files=FALSE)
```

#### Arguments ####

* seuratObject Seurat object following a protocol such as https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
* num_genes Number of genes for the top_n function. Default is 50.
* write_files Save the output files to .txt format. Defauly is FALSE.


#### Value ####

* newExpressionFile - Seurat expression file in ICGS format (ICGS genes)
* newFullExpressionFile - Seurat expression file in ICGS format (all genes)
* newGroupsFile - Groups file ICGS format


<s>
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
</s>

## Seurat and ICGS data: ##

```javascript
Main_Doublet_Decon(rawDataFile, groupsFile, filename, location,
  fullDataFile = NULL, removeCC = FALSE, species = "mmu", rhop = 1,
  write = TRUE, PMF = TRUE, useFull = FALSE, heatmap = TRUE, centroids=FALSE, num_doubs=100, 
  only50=FALSE, min_uniq=4, nCores=-1)
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
* num_doubs: The user defined number of doublets to make for each pair of clusters. Default is 100.
* only50: use only synthetic doublets created with 50%/50% mix of parent cells, as opposed to the extended option of 30%/70% and 70%/30%, default is FALSE.
* min_uniq: minimum number of unique genes required for a cluster to be rescued, default is 4.
* nCores: number of cores to be used during rescue step. Default is -1 for automatically detected.

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

Data for this example can be found in this GitHub repository. Examples are given for both Seurat_Pre_Process() and Improved_Seurat_Pre_Process(), though the latter is prefered if using Seurat 3.

```javascript
location="/Users/xxx/xxx/" #Update as needed 

<s>
#Seurat_Pre_Process()
expressionFile=paste0(location, "counts.txt")
genesFile=paste0(location, "Top50Genes.txt")
clustersFile=paste0(location, "Cluster.txt")
newFiles=Seurat_Pre_Process(expressionFile, genesFile, clustersFile)
</s>

#Improved_Seurat_Pre_Process()
seuratObject=readRDS("seurat.rds")
newFiles=Improved_Seurat_Pre_Process(seuratObject, num_genes=50, write_files=FALSE)

filename="PBMC_example"
write.table(newFiles$newExpressionFile, paste0(location, filename, "_expression"), sep="\t")
write.table(newFiles$newFullExpressionFile, paste0(location, filename, "_fullExpression"), sep="\t")
write.table(newFiles$newGroupsFile, paste0(location, filename , "_groups"), sep="\t", col.names = F)

results=Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                           groupsFile=newFiles$newGroupsFile, 
                           filename=filename, 
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
                           only50=FALSE,
                           min_uniq=4,
                           nCores=-1)
```                           
