#' Main DoubletDecon v1.0.1
#'
#' This is the main function. This function identifies clusters of doublets with a combination of deconvolution analysis and unique gene expression and individual doublet cells with deconvolution analysis.
#' @param rawDataFile Name of file containing ICGS or Seurat expression data (gene by cell)
#' @param groupsFile Name of file containing group assignments (3 column: cell, group(numeric), group(numeric or character))
#' @param filename Unique filename to be incorporated into the names of outputs from the functions.
#' @param location Directory where output should be stored
#' @param fullDataFile Name of file containing full expression data (gene by cell). Default is NULL.
#' @param removeCC Remove cell cycle gene cluster by KEGG enrichment. Default is FALSE.
#' @param species Species as scientific species name, KEGG ID, three letter	species abbreviation, or NCBI ID. Default is "mmu".
#' @param rhop x in mean+x*SD to determine upper cutoff for correlation in the blacklist. Default is 1.
#' @param write Write output files as .txt files. Default is TRUE.
#' @param PMF Use step 2 (unique gene expression) in doublet determination criteria. Default is TRUE.
#' @param useFull Use full gene list for PMF analysis. Requires fullDataFile. Default is FALSE.
#' @param heatmap Boolean value for whether to generate heatmaps. Default is TRUE. Can be slow to datasets larger than ~3000 cells.
#' @param centroids Use centroids as references in deconvolution instead of the default medoids.
#' @param num_doubs The user defined number of doublets to make for each pair of clusters. Default is 30.
#' @param downsample allows for downsampling of cells when using full expression matrix (use with large datasets), default is "none".
#' @param sample_num number of cells per cluster with downsampling with "even", percent of cluster with "prop".
#' @param only50 use only synthetic doublets created with 50%/50% mix of parent cells, as opposed to the extended option of 30%/70% and 70%/30%, default is TRUE.
#' @return data_processed - new expression file (cleaned).
#' @return groups_processed = new groups file (cleaned).
#' @return PMF_results = pseudo marker finder t-test results (gene by cluster).
#' @return DRS_doublet_table = each cell and whether it is called a doublet by deconvolution analysis.
#' @return DRS_results = results of deconvolution analysis (cell by cluster) in percentages.
#' @return Decon_called_freq = percentage of doublets called in each cluster by deconvolution analysis.
#' @return Final_doublets_groups = new groups file containing only doublets.
#' @return Final_nondoublets_groups = new groups file containing only non doublets.
#' @keywords doublet decon main
#' @export

Main_Doublet_Decon<-function(rawDataFile, groupsFile, filename, location, fullDataFile=NULL, removeCC=FALSE, species="mmu", rhop=1, write=TRUE, PMF=TRUE, useFull=FALSE, heatmap=TRUE, centroids=FALSE, num_doubs=30, downsample="none", sample_num=NULL, only50=TRUE){

  #load required packages
  require(DeconRNASeq)
  require(gplots)
  require(dplyr)
  require(MCL)
  require(clusterProfiler)
  require(mygene)
  require(hopach)
  require(as.color)

  #Set up log file
  log_file_name=paste0(location, Sys.time(),".log")
  log_con <- file(log_file_name)
  cat(paste0("filename: ",filename), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("location: ",location), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("removeCC: ",removeCC), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("species: ",species), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("rhop: ",rhop), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("write: ",write), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("PMF: ",PMF), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("useFull: ",useFull), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("heatmap: ",heatmap), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("centroids: ",centroids), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("num_doubs: ",num_doubs), file=log_file_name, append=TRUE, sep="\n")

  #Read in data
  cat("Reading data...", file=log_file_name, append=TRUE, sep="\n")
  cat("Reading data...", sep="\n")
  if(class(rawDataFile)=="character"){
    rawData=read.table(rawDataFile, sep="\t",header=T, row.names=1)
  }else{
    rawData=rawDataFile
  }

  if(class(groupsFile)=="character"){
    groups=read.table(groupsFile, sep="\t",header=F, row.names=1)

  }else{
    groups=groupsFile
  }

  #Clean up data and groups file
  cat("Processing raw data...", file=log_file_name, append=TRUE, sep="\n")
  cat("Processing raw data...", sep="\n")
  data=Clean_Up_Input(rawData, groups, log_file_name=log_file_name)
  og_processed_data=data$processed
  groups=data$groups

  #Centroids or medoids?
  if(centroids==TRUE){
    centroid_flag=TRUE
  }else{
    centroid_flag=FALSE
  }

  #Original data heatmap
  if(heatmap==TRUE){
    cat("Creating original data heatmap...", file=log_file_name, append=TRUE, sep="\n")
    cat("Creating original data heatmap...", sep="\n")
    breaks=seq(0, #start point of color key
               as.numeric(quantile(data.matrix(data$processed[2:nrow(data$processed), 2:ncol(data$processed)]), 0.99)),  #end point of color key
               by=0.05) #length of sub-division
    mycol <- colorpanel(n=length(breaks)-1, low="black", high= "yellow") #heatmap colors
    suppressWarnings(DDheatmap(data.matrix(data$processed[2:nrow(data$processed), 2:ncol(data$processed)]), #the data matrix
                               Colv=FALSE, # No clustering of columns
                               Rowv = FALSE, #no clustering of rows
                               dendrogram = "none", #do not generate dendrogram
                               col=mycol, #colors used in heatmap
                               ColSideColors = as.color(Renumber(data$processed[1,2:ncol(data$processed)]), alpha=1, seed=4), #column color bar
                               RowSideColors = as.color(Renumber(data$processed[2:nrow(data$processed),1]), alpha=1, seed=2), # row color bar
                               breaks=breaks, #color key details
                               trace="none", #no trace on map
                               na.rm=TRUE, #ignore missing values
                               margins = c(5,5), # size and layout of heatmap window
                               labRow=NA, #turn off gene labels
                               labCol=NA, #turn off cell labels
                               xlab = "Samples", #x axis title
                               ylab =  "Genes", # y axis title
                               main = paste0("Original data: ", filename))) #main title
    }

  #Remove cell cycle gene cluster (optional)
  if(removeCC==TRUE){
    cat("Removing cell cycle clusters...", file=log_file_name, append=TRUE, sep="\n")
    cat("Removing cell cycle clusters...", sep="\n")
    data=Remove_Cell_Cycle(data$processed, species, log_file_name)
  }else{
    data=data$processed
  }
  if(write==TRUE){
    write.table(data, paste0(location, "data_processed_", filename, ".txt"), sep="\t")
    write.table(groups, paste0(location, "groups_processed_", filename, ".txt"), sep="\t")
  }

  #Calculate medoids, medoid correlations, blacklist to create new combine medoids
  cat("Combining similar clusters...", file=log_file_name, append=TRUE, sep="\n")
  cat("Combining similar clusters...", sep="\n")
  BL=Blacklist_Groups(data, groups, rhop, centroid_flag, log_file_name)
  newMedoids=BL$newMedoids
  groupsMedoids=BL$newGroups

  #Create synthetic doublets to get average synthetic profiles
  cat("Creating synthetic doublet profiles...", file=log_file_name, append=TRUE, sep="\n")
  cat("Creating synthetic doublet profiles...", sep="\n")
  if(.Platform$OS.type=="unix"){
    sink("/dev/null") #hides DeconRNASeq output
    synthProfilesx=Synthetic_Doublets(data, groups, groupsMedoids, newMedoids, num_doubs, log_file_name=log_file_name)
    sink()
  }else{
    synthProfilesx=Synthetic_Doublets(data, groups, groupsMedoids, newMedoids, num_doubs, log_file_name=log_file_name)
  }
  synthProfiles=synthProfilesx$averagesAverages
  doubletCellsInput2=synthProfilesx$doubletCellsInput2
  if(write==TRUE){
    write.table(doubletCellsInput2, paste0(location, "Synth_doublet_info_", filename, ".txt"), sep="\t")
  }

  #Calculate doublets using DeconRNASeq
  cat("Step 1: Removing possible doublets...", file=log_file_name, append=TRUE, sep="\n")
  cat("Step 1: Removing possible doublets...", sep="\n")
  if(.Platform$OS.type=="unix"){
    sink("/dev/null") #hides DeconRNASeq output
    doubletTable=Is_A_Doublet(data, newMedoids, groups, synthProfiles, log_file_name=log_file_name)
    sink()
  }else{
    doubletTable=Is_A_Doublet(data, newMedoids, groups, synthProfiles, log_file_name=log_file_name)
  }
  if(write==TRUE){
    write.table(doubletTable$isADoublet, paste0(location, "DRS_doublet_table_", filename, ".txt"), sep="\t")
    write.table(doubletTable$resultsreadable, paste0(location, "DRS_results_", filename, ".txt"), sep="\t")
  }

  #Recluster doublets and non-doublets
  cat("Step 2: Re-clustering possible doublets...", file=log_file_name, append=TRUE, sep="\n")
  cat("Step 2: Re-clustering possible doublets...", sep="\n")
  reclusteredData=Recluster(isADoublet=doubletTable$isADoublet, data, groups, log_file_name = log_file_name)
  data=reclusteredData$newData2$processed
  groups=reclusteredData$newData2$groups
  if(write==TRUE){
    write.table(data, paste0(location, "data_processed_reclust_", filename, ".txt"), sep="\t")
    write.table(groups, paste0(location, "groups_processed_reclust_", filename, ".txt"), sep="\t")
  }

  #Run Pseudo Marker Finder to identify clusters with no unique gene expression
  cat("Step 3: Rescuing cells with unique gene expression...", file=log_file_name, append=TRUE, sep="\n")
  cat("Step 3: Rescuing cells with unique gene expression...", sep="\n")
  if(useFull==TRUE){
    if(class(fullDataFile)=="character"){
      full_data=read.table(fullDataFile, sep="\t",header=T, row.names=1)
    }else{
      full_data=fullDataFile
    }
    full_data2=Clean_Up_Input(full_data, groups)$processed
    PMFresults=Pseudo_Marker_Finder(groups, data, full_data2, downsample, sample_num, log_file_name=log_file_name)
  }else{
    PMFresults=Pseudo_Marker_Finder(groups, data, full_data2=NULL, downsample, sample_num, log_file_name=log_file_name)
  }
  if(write==TRUE){
    write.table(PMFresults, paste0(location, "new_PMF_results_", filename, ".txt"), sep="\t")
  }

  #Doublet Detection method 2: Pseudo_Marker_Finder
  allClusters=unique(groups[,1])
  hallmarkClusters=as.numeric(unique(PMFresults[,2]))
  newDoubletClusters=setdiff(allClusters, hallmarkClusters)

  #Doublet Detection method 1: Is_A_Doublet
  uniqueClusters=as.character(unique(groups[,2]))
  DeconCalledFreq=as.data.frame(matrix(nrow=length(allClusters), ncol=1), row.names = uniqueClusters)
  if(recluster!="none"){
    DeconCalledFreq=reclusteredData$decon
  }else{
    for(clus in allClusters){
      temp1=subset(doubletTable$isADoublet, Group_Cluster==uniqueClusters[clus])
      DeconCalledFreq[clus,1]=(length(which(temp1$isADoublet==TRUE))/nrow(temp1))*100
    }
  }

  #Combine to find real doublets
  if(PMF==FALSE){
    finalDoubletClusters=which(DeconCalledFreq>50)
  }else{
    finalDoubletClusters=intersect(which(DeconCalledFreq>50), newDoubletClusters)
  }

  #Results
  finalDoubletCellCall=subset(groups, groups[,1] %in% finalDoubletClusters)
  finalNotDoubletCellCall=subset(groups, !(groups[,1] %in% finalDoubletClusters))
  if(write==TRUE){
    write.table(finalDoubletCellCall, paste0(location, "Final_doublets_groups_", filename, ".txt"), sep="\t")
    write.table(finalNotDoubletCellCall, paste0(location, "Final_nondoublets_groups_", filename, ".txt"), sep="\t")
  }


  #Subset expression matrix for doublets and save
  doublets_matrix=cbind(og_processed_data[,1],og_processed_data[,which(colnames(og_processed_data) %in% row.names(finalDoubletCellCall))])
  if(write==TRUE){
    write.table(doublets_matrix, paste0(location, "Final_doublets_exp_", filename, ".txt"), sep="\t")
  }

  #Heatmap of cells removed as doubets
  if(heatmap==TRUE){
    cat("Creating doublets heatmap...", file=log_file_name, append=TRUE, sep="\n")
    cat("Creating doublets heatmap...", sep="\n")
    breaks=seq(0, #start point of color key
               as.numeric(quantile(data.matrix(doublets_matrix[2:nrow(doublets_matrix), 2:ncol(doublets_matrix)]), 0.99)),  #end point of color key
               by=0.05) #length of sub-division
    mycol <- colorpanel(n=length(breaks)-1, low="black", high= "yellow") #heatmap colors
    suppressWarnings(DDheatmap(data.matrix(doublets_matrix[2:nrow(doublets_matrix), 2:ncol(doublets_matrix)]), #the data matrix
                               Colv=FALSE, # No clustering of columns
                               Rowv = FALSE, #no clustering of rows
                               col=mycol, #colors used in heatmap
                               dendrogram="none", #turn of dendrogram generation
                               ColSideColors = as.color(Renumber(doublets_matrix[1,2:ncol(doublets_matrix)]), alpha=1, seed=4), #column color bar
                               RowSideColors = as.color(Renumber(doublets_matrix[2:nrow(doublets_matrix),1]), alpha=1, seed=2), # row color bar
                               breaks=breaks, #color key details
                               trace="none", #no trace on map
                               na.rm=TRUE, #ignore missing values
                               margins = c(5,5), # size and layout of heatmap window
                               labRow=NA, #turn off gene labels
                               labCol=NA, #turn off cell labels
                               xlab = "Samples", #x axis title
                               ylab =  "Genes", # y axis title
                               main = paste0("Doublets: ", filename))) #main title)
  }

  #Subset expression matrix for non-doublets and save
  nondoublets_matrix=cbind(og_processed_data[,1],og_processed_data[,which(colnames(og_processed_data) %in% row.names(finalNotDoubletCellCall))])
  if(write==TRUE){
    write.table(nondoublets_matrix, paste0(location, "Final_nondoublets_exp_", filename, ".txt"), sep="\t")
  }

  #New heatmap of non-doublet cells
  if(heatmap==TRUE){
    cat("Creating non-doublets heatmap...", file=log_file_name, append=TRUE, sep="\n")
    cat("Creating non-doublets heatmap...", sep="\n")
    breaks=seq(0, #start point of color key
               as.numeric(quantile(data.matrix(nondoublets_matrix[2:nrow(nondoublets_matrix), 2:ncol(nondoublets_matrix)]), 0.99)),  #end point of color key
               by=0.05) #length of sub-division
    mycol <- colorpanel(n=length(breaks)-1, low="black", high= "yellow") #heatmap colors
    suppressWarnings(DDheatmap(data.matrix(nondoublets_matrix[2:nrow(nondoublets_matrix), 2:ncol(nondoublets_matrix)]), #the data matrix
                               Colv=FALSE, # No clustering of columns
                               Rowv = FALSE, #no clustering of rows
                               col=mycol, #colors used in heatmap
                               dendrogram="none", #turn of dendrogram generation
                               ColSideColors = as.color(Renumber(nondoublets_matrix[1,2:ncol(nondoublets_matrix)]), alpha=1, seed=4), #column color bar
                               RowSideColors = as.color(Renumber(nondoublets_matrix[2:nrow(nondoublets_matrix),1]), alpha=1, seed=2), # row color bar
                               breaks=breaks, #color key details
                               trace="none", #no trace on map
                               na.rm=TRUE, #ignore missing values
                               margins = c(5,5), # size and layout of heatmap window
                               labRow=NA, #turn off gene labels
                               labCol=NA, #turn off cell labels
                               xlab = "Samples", #x axis title
                               ylab =  "Genes", # y axis title
                               main = paste0("Non-Doublets: ", filename))) #main title
  }

  #close the log file connection
  close(log_con)

  return(list(data_processed=data,
              groups_processed=groups,
              DRS_doublet_table=doubletTable$isADoublet,
              DRS_results=doubletTable$resultsreadable,
              PMF_results=PMFresults,
              Final_doublets_groups=finalDoubletCellCall,
              Final_nondoublets_groups=finalNotDoubletCellCall,
              Final_doublets_exp=doublets_matrix,
              Final_nondoublets_exp=nondoublets_matrix,
              Synth_doublet_info=doubletCellsInput2))

}
