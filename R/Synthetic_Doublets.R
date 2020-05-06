#' Synthetic Doublets
#'
#' This function creates synthetic doublets by averaging gene expression profiles from each combination of clusters to generate deconvolution profiles for each type of doublet.
#' @param data Processed data from Clean_Up_Input (or Remove_Cell_Cycle).
#' @param groups Processed groups file from Clean_Up_Input.
#' @param groupsMedoids New groups file based on blacklisted clusters from Blacklist_Groups
#' @param newMedoids New combined medoids from Blacklist_groups.
#' @param num_doubs The user defined number of doublets to make for each pair of clusters
#' @param log_file_name used for saving run notes to log file
#' @param only50 use only synthetic doublets created with 50\%/50\% mix of parent cells, as opposed to the extended option of 30\%/70\% and 70\%/30\%, default is TRUE.
#' @return averagesAverages = average deconvolution profiles for each combination of cell types.
#' @return doubletCellsInput2 = deconvolution profiles for synthetic doublet cells (for quality control).
#' @keywords synthetic
#' @export

Synthetic_Doublets<-function(data, groups, groupsMedoids, newMedoids, num_doubs, log_file_name, only50, location){

  #Override the original groups for making synthetics with groups based on the blacklisted clusters
  groups=groupsMedoids

  ## NEW CODE 28Feb2020 ##
  prop_doubs=num_doubs #10 percent of the total number of cells
  ndub=round(ncol(data)*prop_doubs)
  doubletCellsInput=as.data.frame(matrix(ncol=2, nrow=(ndub)))
  doubletCellsInput2=as.data.frame(matrix(ncol=4, nrow=(ndub))) #for info
  for(pair in 1:ndub){
    sampled_pair=sample(row.names(groups), 2, replace=FALSE)
    sampled_pair_groups=groups[which(row.names(groups) %in% sampled_pair),2]
    while(sampled_pair_groups[1]==sampled_pair_groups[2]){ #make sure that they aren't sampled from the same group
      sampled_pair=sample(row.names(groups), 2, replace=FALSE)
      sampled_pair_groups=groups[which(row.names(groups) %in% sampled_pair),2]
    }
    doubletCellsInput[pair, 1]=sampled_pair[1]
    doubletCellsInput[pair, 2]=sampled_pair[2]
    doubletCellsInput2[pair, 1]=sampled_pair[1]
    doubletCellsInput2[pair, 2]=sampled_pair[2]
    doubletCellsInput2[pair, 3]=sampled_pair_groups[1]
    doubletCellsInput2[pair, 4]=sampled_pair_groups[2]
  }
  
  ## SECTION REPLACED 28Feb2020 ##
  # #Number of doublets created per cluster pair
  # ndub=num_doubs
  # 
  # #Make data frame to hold new doublets (30 per combination of clusters)
  # pairs=combn(unique(groups[,2]), 2)
  # doubletCellsInput=as.data.frame(matrix(ncol=2, nrow=((length(pairs)/2)*ndub)))
  # doubletCellsInput2=as.data.frame(matrix(ncol=4, nrow=((length(pairs)/2)*ndub))) #for info
  # 
  # #For each combination of medoids, select 30 pairs
  # for(pair in 0:((length(pairs)/2)-1)){
  #   for(synth in 1:ndub){
  #     doubletCellsInput[synth+(ndub*pair),1]=sample(row.names(subset(groups, groups[,2]==pairs[1,(pair+1)])),1,replace=FALSE)
  #     doubletCellsInput2[synth+(num_doubs*pair),1]=doubletCellsInput[synth+(num_doubs*pair),1]
  #     doubletCellsInput[synth+(ndub*pair),2]=sample(row.names(subset(groups, groups[,2]==pairs[2,(pair+1)])),1,replace=FALSE)
  #     doubletCellsInput2[synth+(num_doubs*pair),2]=doubletCellsInput[synth+(num_doubs*pair),2]
  #     doubletCellsInput2[synth+(num_doubs*pair),3]=pairs[1,(pair+1)]
  #     doubletCellsInput2[synth+(num_doubs*pair),4]=pairs[2,(pair+1)]
  #   }
  # }

  doubletAverages=as.data.frame(matrix(ncol=nrow(doubletCellsInput), nrow=nrow(data)))
  row.names(doubletAverages)=row.names(data)
  doubletAverages[1,]=rep((length(unique(groups[,1]))+1), ncol(doubletAverages))

  doubletAverages=as.data.frame(matrix(ncol=nrow(doubletCellsInput)*3, nrow=nrow(data)))
  row.names(doubletAverages)=row.names(data)

  doubletAverages[1,]=rep((length(unique(groups[,1]))+1), ncol(doubletAverages)*3)

  for(doublet in 1:(ncol(doubletAverages)/3)){
    cell1=as.character(doubletCellsInput[doublet,1])
    cell2=as.character(doubletCellsInput[doublet,2])
    expression1=data[2:nrow(data),which(colnames(data)==cell1)]
    expression2=data[2:nrow(data),which(colnames(data)==cell2)]
    temp=cbind(expression1, expression2)
    newExpression=apply(temp, 1, weighted.mean, c(0.5,0.5))
    newExpression_a=apply(temp, 1, weighted.mean, c(0.70, 0.30))
    newExpression_b=apply(temp, 1, weighted.mean, c(0.30, 0.70))
    doubletAverages[2:nrow(doubletAverages), doublet]=newExpression
    doubletAverages[2:nrow(doubletAverages), doublet+(ncol(doubletAverages)/3)]=newExpression_a
    doubletAverages[2:nrow(doubletAverages), doublet+((ncol(doubletAverages)/3)*2)]=newExpression_b
    colnames(doubletAverages)[doublet]=paste0(cell1,"-", cell2, "-even")
    colnames(doubletAverages)[doublet+(ncol(doubletAverages)/3)]=paste0(cell1,"-", cell2, "-one")
    colnames(doubletAverages)[doublet+((ncol(doubletAverages)/3)*2)]=paste0(cell1,"-", cell2, "-two")
  }

  #50%/50% references only or 30%/70% and 70%/30% included
  if(only50==TRUE){
    mult=1
  }else{
    mult=3
  }

  write.table(doubletAverages, paste0("~/Documents/Projects/DoubletDetection/DF_DD/", "exp_synth_results_testing_T.txt"), sep="\t")
  write.table(data[-1,-1], paste0("~/Documents/Projects/DoubletDetection/DF_DD/", "exp_real_results_testing_T.txt"), sep="\t")
  expAverages=data[-1,-1]
  results=DeconRNASeq(doubletAverages[2:nrow(doubletAverages),], newMedoids)
  resultsreadable=round(results$out.all*100,2)
  write.table(resultsreadable, paste0(location, "resultsreadable_synths.txt"), sep="\t")
  write.table(resultsreadable, paste0("~/Documents/Projects/DoubletDetection/DF_DD/", "DRS_synth_results_testing.txt"), sep="\t")
  row.names(resultsreadable)=colnames(doubletAverages)

  ## CODE HIDDEN 28Feb2020 ##
  # averagesAverages=as.data.frame(matrix(ncol=ncol(resultsreadable), nrow=(length(pairs)/2)*mult))
  # colnames(averagesAverages)=colnames(resultsreadable)
  # i=1
  # for(clust in 1:nrow(averagesAverages)){
  #   averagesAverages[clust,]=apply(resultsreadable[i:(i+(num_doubs-1)),], 2, mean)
  #   if(clust %in% 1:(length(pairs)/2)){
  #     row.names(averagesAverages)[clust]=paste0(pairs[1,clust], "-", pairs[2,clust], "-even")
  #   }else if(clust %in% ((length(pairs)/2)+1):((length(pairs)/2)*2)){
  #     row.names(averagesAverages)[clust]=paste0(pairs[1,clust-(length(pairs)/2)], "-", pairs[2,clust-(length(pairs)/2)], "-one")
  #   }else{
  #     row.names(averagesAverages)[clust]=paste0(pairs[1,clust-length(pairs)], "-", pairs[2,clust-length(pairs)], "-two")
  #   }
  # 
  #   i=i+num_doubs
  # }
  averagesAverages=NULL

  #return(list(averagesAverages=averagesAverages, doubletCellsInput2=doubletCellsInput2))
  return(list(averagesAverages=averagesAverages, doubletCellsInput2=doubletCellsInput2, synth_DRS_results=resultsreadable))
}

