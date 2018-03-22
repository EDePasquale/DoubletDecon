#' Synthetic Doublets
#'
#' This function creates synthetic doublets by averaging gene expression profiles from each combination of clusters to generate deconvolution profiles for each type of doublet.
#' @param data Processed data from Clean_Up_Input (or Remove_Cell_Cycle).
#' @param groups Processed groups file from Clean_Up_Input.
#' @param groupsMedoids New groups file based on blacklisted clusters from Blacklist_Groups
#' @param newMedoids New combined medoids from Blacklist_groups.
#' @return averagesAverages - average deconvolution profiles for each combination of cell types.
#' @keywords synthetic
#' @export

Synthetic_Doublets<-function(data, groups, groupsMedoids, newMedoids){

  #Override the original groups for making synthetics with groups based on the blacklisted clusters
  groups=groupsMedoids

  ndub=30
  #Make data frame to hold new doublets (30 per combination of clusters)
  pairs=combn(unique(groups[,2]), 2)
  doubletCellsInput=as.data.frame(matrix(ncol=2, nrow=((length(pairs)/2)*ndub)))

  #For each combination of medoids, select 30 pairs

  for(pair in 0:((length(pairs)/2)-1)){
    for(synth in 1:ndub){
      doubletCellsInput[synth+(ndub*pair),1]=sample(row.names(subset(groups, groups[,2]==pairs[1,(pair+1)])),1,replace=FALSE)
      doubletCellsInput[synth+(ndub*pair),2]=sample(row.names(subset(groups, groups[,2]==pairs[2,(pair+1)])),1,replace=FALSE)

    }
  }

  doubletAverages=as.data.frame(matrix(ncol=nrow(doubletCellsInput), nrow=nrow(data)))
  row.names(doubletAverages)=row.names(data)
  doubletAverages[1,]=rep((length(unique(groups[,1]))+1), ncol(doubletAverages))

  for(doublet in 1:ncol(doubletAverages)){
    cell1=as.character(doubletCellsInput[doublet,1])
    cell2=as.character(doubletCellsInput[doublet,2])

    expression1=data[2:nrow(data),which(colnames(data)==cell1)]
    expression2=data[2:nrow(data),which(colnames(data)==cell2)]
    temp=cbind(expression1, expression2)

    newExpression=apply(temp, 1, mean)

    doubletAverages[2:nrow(doubletAverages), doublet]=newExpression

    colnames(doubletAverages)[doublet]=paste0(cell1,"-", cell2)
  }

  #write.table(doubletAverages, paste0(location, "new_synths.txt"), sep="\t")
  results=DeconRNASeq(doubletAverages[2:nrow(doubletAverages),], newMedoids)
  resultsreadable=round(results$out.all*100,2)
  row.names(resultsreadable)=colnames(doubletAverages)

  averagesAverages=as.data.frame(matrix(ncol=ncol(resultsreadable), nrow=length(pairs)/2))
  colnames(averagesAverages)=colnames(resultsreadable)
  kitty=1
  for(clust in 1:(length(pairs)/2)){
    averagesAverages[clust,]=apply(resultsreadable[kitty:(kitty+(ndub-1)),], 2, mean)
    row.names(averagesAverages)[clust]=paste0(pairs[1,clust], "-", pairs[2,clust])
    kitty=kitty+ndub
  }
  return(averagesAverages)
}
