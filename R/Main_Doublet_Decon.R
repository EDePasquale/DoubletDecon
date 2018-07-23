#' Main DoubletDecon
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
#' @param recluster Recluster deconvolution classified doublets and non-doublets seperately using hopach or deconvolution classifications.
#' @param PMF Use step 2 (unique gene expression) in doublet determination criteria. Default is TRUE.
#' @param useFull Use full gene list for PMF analysis. Requires fullDataFile. Default is FALSE.
#' @param heatmap Boolean value for whether to generate heatmaps. Default is TRUE. Can be slow to datasets larger than ~3000 cells.
#' @param centroids Use centroids as references in deconvolution instead of the default medoids.
#' @param num_doubs The user defined number of doublets to make for each pair of clusters. Default is 30.
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

#Main_Doublet_Decon<-function(rawDataFile, groupsFile, filename, removeCC=FALSE, species="mmu", rhop=1, write=TRUE, recluster="none", isADoublet=1, PMF=TRUE){
Main_Doublet_Decon<-function(rawDataFile, groupsFile, filename, location, fullDataFile=NULL, removeCC=FALSE, species="mmu", rhop=1, write=TRUE, recluster="doublets_decon", PMF=TRUE, useFull=FALSE, heatmap=TRUE, centroids=FALSE, num_doubs=30){

  #load required packages
  require(DeconRNASeq)
  require(gplots)
  require(dplyr)
  require(MCL)
  require(clusterProfiler)
  require(mygene)
  require(hopach)
  require(as.color)

  #Read in data
  print("Reading data...")
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
  print("Processing raw data...")
  data=Clean_Up_Input(rawData, groups)
  og_processed_data=data$processed
  rowgroups=data$groups

  #TODO: work in progress
  #Check for sparsity
  #print("Checking for sparsity...")
  # centroid_flag=FALSE
  # if(!is.na(data$processed[2,1])){
  #   temp=apply(data$processed[2:nrow(data$processed), 2:ncol(data$processed)], 1, median)
  #   temp2=cbind(temp, data$processed[2:nrow(data$processed), 1])
  #   for(geneclust in 1:length(unique(temp2[,2]))){ #by cluster, if a gene cluster has a high amount of 0s for its mediod then move all medoids to centroids
  #     temp3=temp2[which(temp2[,2]==geneclust),]
  #     temp4=length(which(temp3[,1]==0))/length(temp3[,1])
  #     if(temp4>0.9){
  #       centroid_flag=TRUE
  #     }
  #   }
  # }
  # if(centroid_flag==TRUE){
  #   print("High sparsity in dataset, moving to centroids...")
  # }

  if(centroids==TRUE){
    centroid_flag=TRUE
  }else{
    centroid_flag==FALSE
  }

  #Original data heatmap
  if(heatmap==TRUE){
    print("Creating original data heatmap...")
    breaks=seq(0, #start point of color key
               as.numeric(quantile(data.matrix(data$processed[2:nrow(data$processed), 2:ncol(data$processed)]), 0.99)),  #end point of color key
               by=0.05) #length of sub-division
    mycol <- colorpanel(n=length(breaks)-1, low="black", high= "yellow") #heatmap colors
    suppressWarnings(heatmap.2(data.matrix(data$processed[2:nrow(data$processed), 2:ncol(data$processed)]), #the data matrix
                               Colv=FALSE, # No clustering of columns
                               Rowv = FALSE, #no clustering of rows
                               col=mycol, #colors used in heatmap
                               ColSideColors = as.color(Renumber(data$processed[1,2:ncol(data$processed)]), alpha=1, seed=4), #column color bar
                               RowSideColors = as.color(Renumber(data$processed[2:nrow(data$processed),1]), alpha=1, seed=2), # row color bar
                               breaks=breaks, #color key details
                               trace="none", #no trace on map
                               na.rm=TRUE, #ignore missing values
                               margins = c(5,5), # size and layout of heatmap window
                               xlab = "Samples", #x axis title
                               ylab =  "Genes", # y axis title
                               main = paste0("Original data: ", filename))) #main title
  }

  #Remove cell cycle gene cluster (optional)
  if(removeCC==TRUE){
    print("Removing cell cycle clusters...")
    data=Remove_Cell_Cycle(data$processed, species)
  }else{
    data=data$processed
  }

  if(write==TRUE){
    write.table(data, paste0(location, "data_processed_", filename, ".txt"), sep="\t")
    write.table(groups, paste0(location, "groups_processed_", filename, ".txt"), sep="\t")
  }

  #Calculate medoids, medoid correlations, blacklist to create new combine medoids
  print("Combining similar clusters...")
  BL=Blacklist_Groups(data, groups, rhop, centroid_flag)
  newMedoids=BL$newMedoids
  groupsMedoids=BL$newGroups

  #Create synthetic doublets to get average synthetic profiles
  print("Creating synthetic doublet profiles...")
  sink("/dev/null") #hides DeconRNASeq output
  synthProfiles=Synthetic_Doublets(data, groups, groupsMedoids, newMedoids, num_doubs)
  sink()

  #Calculate doublets using DeconRNASeq
  print("Step 1: Removing possible doublets...")
  sink("/dev/null") #hides DeconRNASeq output
  doubletTable=Is_A_Doublet(data, newMedoids, groups, synthProfiles)
  sink()
  if(write==TRUE){
    write.table(doubletTable$isADoublet, paste0(location, "DRS_doublet_table_", filename, ".txt"), sep="\t")
    write.table(doubletTable$resultsreadable, paste0(location, "DRS_results_", filename, ".txt"), sep="\t")
  }

  #Recluster doublets and non-doublets
  print("Step 2: Re-clustering possible doublets...")
  reclusteredData=Recluster(doubletTable$isADoublet, data, recluster, groups)
  data=reclusteredData$newData2$processed
  groups=reclusteredData$newData2$groups
    if(write==TRUE){
      write.table(data, paste0(location, "data_processed_reclust_", filename, ".txt"), sep="\t")
      write.table(groups, paste0(location, "groups_processed_reclust_", filename, ".txt"), sep="\t")
    }


  #Run Pseudo Marker Finder to identify clusters with no unique gene expression
  print("Step 3: Rescuing cells with unique gene expression...")
  if(useFull==TRUE){
    if(class(fullDataFile)=="character"){
      full_data=read.table(fullDataFile, sep="\t",header=T, row.names=1)
    }else{
      full_data=fullDataFile
    }
    full_data2=Clean_Up_Input(full_data, groups)$processed
    PMFresults=Pseudo_Marker_Finder(groups, data, full_data2) #TODO: Not sure what to do about the groups here...
  }else{
    PMFresults=Pseudo_Marker_Finder(groups, data, full_data2=NULL)
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


  if(write==TRUE){
    write.table(DeconCalledFreq, paste0(location, "Decon_called_freq_", filename, ".txt"), sep="\t")
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
    print("Creating doublets heatmap...")
    breaks=seq(0, #start point of color key
               as.numeric(quantile(data.matrix(doublets_matrix[2:nrow(doublets_matrix), 2:ncol(doublets_matrix)]), 0.99)),  #end point of color key
               by=0.05) #length of sub-division
    mycol <- colorpanel(n=length(breaks)-1, low="black", high= "yellow") #heatmap colors
    suppressWarnings(heatmap.2(data.matrix(doublets_matrix[2:nrow(doublets_matrix), 2:ncol(doublets_matrix)]), #the data matrix
                               Colv=FALSE, # No clustering of columns
                               Rowv = FALSE, #no clustering of rows
                               col=mycol, #colors used in heatmap
                               ColSideColors = as.color(Renumber(doublets_matrix[1,2:ncol(doublets_matrix)]), alpha=1, seed=4), #column color bar
                               RowSideColors = as.color(Renumber(doublets_matrix[2:nrow(doublets_matrix),1]), alpha=1, seed=2), # row color bar
                               breaks=breaks, #color key details
                               trace="none", #no trace on map
                               na.rm=TRUE, #ignore missing values
                               margins = c(5,5), # size and layout of heatmap window
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
    print("Creating non-doublets heatmap...")
    breaks=seq(0, #start point of color key
               as.numeric(quantile(data.matrix(nondoublets_matrix[2:nrow(nondoublets_matrix), 2:ncol(nondoublets_matrix)]), 0.99)),  #end point of color key
               by=0.05) #length of sub-division
    mycol <- colorpanel(n=length(breaks)-1, low="black", high= "yellow") #heatmap colors
    suppressWarnings(heatmap.2(data.matrix(nondoublets_matrix[2:nrow(nondoublets_matrix), 2:ncol(nondoublets_matrix)]), #the data matrix
                               Colv=FALSE, # No clustering of columns
                               Rowv = FALSE, #no clustering of rows
                               col=mycol, #colors used in heatmap
                               ColSideColors = as.color(Renumber(nondoublets_matrix[1,2:ncol(nondoublets_matrix)]), alpha=1, seed=4), #column color bar
                               RowSideColors = as.color(Renumber(nondoublets_matrix[2:nrow(nondoublets_matrix),1]), alpha=1, seed=2), # row color bar
                               breaks=breaks, #color key details
                               trace="none", #no trace on map
                               na.rm=TRUE, #ignore missing values
                               margins = c(5,5), # size and layout of heatmap window
                               xlab = "Samples", #x axis title
                               ylab =  "Genes", # y axis title
                               main = paste0("Non-Doublets: ", filename))) #main title
  }

  return(list(data_processed=data,
              groups_processed=groups,
              DRS_doublet_table=doubletTable$isADoublet,
              DRS_results=doubletTable$resultsreadable,
              PMF_results=PMFresults,
              Decon_called_freq=DeconCalledFreq,
              Final_doublets_groups=finalDoubletCellCall,
              Final_nondoublets_groups=finalNotDoubletCellCall,
              Final_doublets_exp=doublets_matrix,
              Final_nondoublets_exp=nondoublets_matrix))

}
