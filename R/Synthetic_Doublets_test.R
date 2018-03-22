#' Synthetic Doublets (for testing)
#'
#' This function creates synthetic doublets by averaging gene expression profiles from each combination of clusters to generate deconvolution profiles for each type of doublet. Doublets should be incorporated into data with cellHarmony prior to testing.
#' @param rawDataFile Name of file containing expression data (gene by cell)
#' @param groupsFile Name of file containing group assignments (3 column: cell, group(numeric), group(numeric or character))
#' @param filename Unique filename to be incorporated into the names of outputs from the functions
#' @param removeCC Remove cell cycle gene cluster by KEGG enrichment. Default is FALSE.
#' @param species Species as scientific species name, KEGG ID, three letter	species abbreviation, or NCBI ID. Default is "mmu".
#' @param corrCutoff x in mean*x*SD to determine upper cutoff for correlation in the blacklist. Default is 1.
#' @param numDubPercent number of cells to be created as doublets based on percentage of total cells. Default is 0.15.
#' @keywords synthetic
#' @export

Synthetic_Doublets_test<-function(rawDataFile, groupsFile, filename, removeCC=FALSE, species="mmu", corrCutoff=1, numDubPercent=0.15){
  
  #load required packages
  library(DeconRNASeq)
  library(gplots)
  library(dplyr)
  library(MCL)
  library(clusterProfiler)
  library(mygene)
  library(multtest)
  library(hopach)
  library(as.color)
  
  #Read in data
  rawData=read.table(rawDataFile, sep="\t",header=T, row.names=1)
  groups=read.table(groupsFile, sep="\t",header=F, row.names=1)
  
  #Clean up data and groups file
  data=Clean_Up_Input(rawData, groups)
  groups=data$groups
  
  #Remove cell cycle gene cluster (optional)
  if(removeCC==TRUE){
    data=Remove_Cell_Cycle(data$processed, species)
  }else{
    data=data$processed
  }
  
  #Calculate medoids, medoid correlations, blacklist to create new combine medoids
  BL=Blacklist_Groups(data, groups, corrCutoff)
  newMedoids=BL$newMedoids
  groupsMedoids=BL$newGroups
  
  #Create data.frame to hold synthetic doublets
  stripRaw=data[2:nrow(data), 2:ncol(data)]
  numSynth=round(ncol(stripRaw) * numDubPercent)
  uniqueBLMedoids=unique(groupsMedoids[,2])
  percentCells=table(groupsMedoids[,2])/nrow(groupsMedoids)*100
  synthRawData=data.frame(matrix(ncol=1, nrow=nrow(stripRaw)), row.names = row.names(stripRaw))
  synthRawData[,1]=row.names(synthRawData)
  
  #Create data.frame to hold synthetic doublet origin information
  testDF=as.data.frame(matrix(ncol=4, nrow=numSynth))
  
  #Create doublets
  for(cell in 1:numSynth){
    
    #choose 2 BL clusters
    twoMedoids=sample(sort(uniqueBLMedoids), 2, replace=FALSE, prob=as.numeric(percentCells))
    testDF[cell,1:2]=twoMedoids
    
    #choose 2 cells
    cell1=sample(row.names(subset(groupsMedoids, groupsMedoids$X3==twoMedoids[1])),1,replace=FALSE)
    cell2=sample(row.names(subset(groupsMedoids, groupsMedoids$X3==twoMedoids[2])),1,replace=FALSE)
    testDF[cell,3]=cell1
    testDF[cell,4]=cell2
    newframe=data.frame(matrix(nrow=nrow(stripRaw), ncol=2), row.names=row.names(stripRaw))
    newframe[,1]=stripRaw[,which(colnames(stripRaw)==cell1)]
    newframe[,2]=stripRaw[,which(colnames(stripRaw)==cell2)]
    
    #sum and average
    sumOfCells=apply(newframe, 1, sum)
    newframe=cbind(newframe,sumOfCells)
    synthRawData[,cell]=newframe[,3]/2
  }

  #Add names and attach to original data
  colnames(synthRawData)=paste0("synth", 1:(numSynth))
  newRawData=cbind(stripRaw, synthRawData)
  
  #Make new groups file (real cells + synthetic doublets)
  groupsSynth=as.data.frame(matrix(nrow=numSynth, ncol=2))
  row.names(groupsSynth)=colnames(synthRawData)
  groupsSynth[,1]=rep(length(unique(groups[,2]))+1, nrow(groupsSynth))
  groupsSynth[,2]=rep("Synth", nrow(groupsSynth))
  colnames(groupsSynth)=colnames(groups)
  groupsSynth=rbind(groups, groupsSynth)
  
  #Add gene clusters
  kitty=rbind(rep(length(unique(groups[,2]))+1, nrow(groupsSynth)), synthRawData)
  kitty2=cbind(data, kitty)
  row.names(kitty2)[1]="column_clusters-flat"
  colnames(kitty2)[1]="row_clusters.flat"
  

  #Write output
  write.table(data.frame("uid"=rownames(synthRawData),synthRawData), paste0(location, "synthExp.", filename, ".txt"), sep="\t", row.names=F)
  #write.table(kitty2, paste0(location, "data_Synth_", filename, ".txt"), sep="\t")
  #write.table(groupsSynth, paste0(location, "groups_Synth_", filename, ".txt"), sep="\t", col.names=FALSE)
  write.table(testDF, paste0(location, "cellsUsed_Synth_", filename, ".txt"), sep="\t")
  
}
