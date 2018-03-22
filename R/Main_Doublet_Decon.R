#' Main DoubletDecon
#'
#' This is the main function. This function identifies clusters of doublets with a combination of deconvolution analysis and unique gene expression and individual doublet cells with deconvolution analysis.
#' @param rawDataFile Name of file containing expression data (gene by cell)
#' @param groupsFile Name of file containing group assignments (3 column: cell, group(numeric), group(numeric or character))
#' @param filename Unique filename to be incorporated into the names of outputs from the functions
#' @param removeCC Remove cell cycle gene cluster by KEGG enrichment. Default is FALSE.
#' @param species Species as scientific species name, KEGG ID, three letter	species abbreviation, or NCBI ID. Default is "mmu".
#' @param corrCutoff x in mean*x*SD to determine upper cutoff for correlation in the blacklist. Default is 1.
#' @param write Write output files as .txt files. Default is TRUE.
#' @param recluster Recluster deconvolution classified doublets and non-doublets deperately using hopach. Default is FALSE.
#' @param isADoublet Use Decon Correlation method (1) or Two Percentages method (2). Default is 2.
#' @return data_processed - new expression file (cleaned).
#' @return groups_processed = new groups file (cleaned).
#' @return PMF_results_1 = pseudo marker finder t statistics (gene by cluster).
#' @return PMF_results_2 = pseudo marker finder gene assignments (which cluster has the highest t-stat for each gene).
#' @return DRS_doublet_table = each cell and whether it is called a doublet by deconvolution analysis.
#' @return DRS_results = results of deconvolution analysis (cell by cluster) in percentages.
#' @return Decon_called_freq = percentage of doublets called in each cluster by deconvolution analysis.
#' @return Final_doublets_groups = new groups file containing only doublets.
#' @return Final_nondoublets_groups = new groups file containing only non doublets.
#' @keywords doublet decon main
#' @export

Main_Doublet_Decon<-function(rawDataFile, groupsFile, filename, removeCC=FALSE, species="mmu", corrCutoff=1, write=TRUE, recluster="none", isADoublet=1){

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

  if(write==TRUE){
    write.table(data, paste0(location, "data_processed_", filename, ".txt"), sep="\t")
    write.table(groups, paste0(location, "groups_processed_", filename, ".txt"), sep="\t")
  }

  #Calculate medoids, medoid correlations, blacklist to create new combine medoids
  BL=Blacklist_Groups(data, groups, corrCutoff)
  newMedoids=BL$newMedoids
  groupsMedoids=BL$newGroups

  #Create synthetic doublets to get average synthetic profiles
  synthProfiles=Synthetic_Doublets(data, groups, groupsMedoids, newMedoids)

  #Calculate doublets using DeconRNASeq
  if(isADoublet==1){
    doubletTable=Is_A_Doublet(data, newMedoids, groups, synthProfiles)
  }else{
    doubletTable=Is_A_Doublet2(data, newMedoids, groups, synthProfiles)
  }
  if(write==TRUE){
    write.table(doubletTable$isADoublet, paste0(location, "DRS_doublet_table_", filename, ".txt"), sep="\t")
    write.table(doubletTable$resultsreadable, paste0(location, "DRS_results_", filename, ".txt"), sep="\t")
  }

  #Recluster doublets and non-doublets with HOPACH (optional)
  if(recluster!="none"){
    print(recluster)
    reclusteredData=Recluster(doubletTable$isADoublet, data, recluster, groups)
    if(length(reclusteredData)<2){
      recluster=reclusteredData
      print("in the if")
      print(recluster)
      #do nothing
    }else{
      print("in the else")
      data=reclusteredData$newData2$processed
      groups=reclusteredData$newData2$groups
      if(write==TRUE){
        write.table(data, paste0(location, "data_processed_reclust_", filename, ".txt"), sep="\t")
        write.table(groups, paste0(location, "groups_processed_reclust_", filename, ".txt"), sep="\t")
      }
    }
  }


  #Run Psuedo Marker Finder to identify clusters with no unique gene expression
  PMFresults=Pseudo_Marker_Finder(groups, data)
  if(write==TRUE){
    write.table(PMFresults$hallmarkTable, paste0(location, "PMF_results_1_", filename, ".txt"), sep="\t")
    write.table(PMFresults$hallmarkTable2, paste0(location, "PMF_results_2_", filename, ".txt"), sep="\t")
  }

  #Doublet Detection method 1: Pseudo_Marker_Finder
  allClusters=1:ncol(PMFresults$hallmarkTable)
  hallmarkClusters=unique(PMFresults$hallmarkTable2[,1])
  newDoubletClusters=setdiff(allClusters, hallmarkClusters)

  #Doublet Detection method 2: Is_A_Doublet
  uniqueClusters=as.character(unique(groups[,2]))
  DeconCalledFreq=as.data.frame(matrix(nrow=length(allClusters), ncol=1), row.names = uniqueClusters)
  if(recluster!="none"){
    print("in the later")
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
  finalDoubletClusters=intersect(which(DeconCalledFreq>50), newDoubletClusters)

  #Results
  if(length(finalDoubletClusters>0)){
    print(paste0(length(finalDoubletClusters), " doublet cluster(s) identified:"))
    print(finalDoubletClusters)
  }else{
    print("0 doublet clusters identified")
  }
  finalDoubletCellCall=subset(groups, groups[,1] %in% finalDoubletClusters)
  finalNotDoubletCellCall=subset(groups, !(groups[,1] %in% finalDoubletClusters))
  if(write==TRUE){
    write.table(finalDoubletCellCall, paste0(location, "Final_doublets_groups_", filename, ".txt"), sep="\t")
    write.table(finalNotDoubletCellCall, paste0(location, "Final_nondoublets_groups_", filename, ".txt"), sep="\t")
  }


  return(list(data_processed=data,
              groups_processed=groups,
              DRS_doublet_table=doubletTable$isADoublet,
              DRS_results=doubletTable$resultsreadable,
              PMF_results_1=PMFresults$hallmarkTable,
              PMF_results_2=PMFresults$hallmarkTable2,
              Decon_called_freq=DeconCalledFreq,
              Final_doublets_groups=finalDoubletCellCall,
              Final_nondoublets_groups=finalNotDoubletCellCall))

}
