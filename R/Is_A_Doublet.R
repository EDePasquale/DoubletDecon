#' Is A Doublet
#'
#' This function uses deconvolution analysis (DeconRNASeq) to evaluate each cell for equal contribution from blacklisted clusters.
#' @param data Processed data from Clean_Up_Input (or Remove_Cell_Cycle).
#' @param newMedoids New combined medoids from Blacklist_groups.
#' @param groups Processed groups file from Clean_Up_Input.
#' @param synthProfiles Average profiles of synthetic doublets from Synthetic_Doublets.
#' @return isADoublet - data.frame with each cell as a row and whether it is called a doublet by deconvolution analysis.
#' @return resultsreadable - data.frame with results of deconvolution analysis (cell by cluster) in percentages.
#' @keywords doublet deconvolution decon
#' @export

Is_A_Doublet<-function(data, newMedoids, groups, synthProfiles){

  #create data frame to store doublets table
  isADoublet=data.frame(matrix(ncol=4,nrow=(ncol(data)-1)))
  rownames(isADoublet)=colnames(data)[2:ncol(data)]
  rownames(newMedoids)=rownames(data)[2:nrow(data)]

  #run DeconRNASeq with new medoids and data
  results=DeconRNASeq(data[2:nrow(data), 2:ncol(data)], newMedoids)
  resultsreadable=round(results$out.all*100,2)
  rownames(resultsreadable)=rownames(isADoublet) #make an easily readable results table

  #get average profiles for cell clusters
  averagesReal=as.data.frame(matrix(ncol=ncol(resultsreadable), nrow=length(unique(groups[,2]))))
  colnames(averagesReal)=colnames(resultsreadable)
  for(clust in 1:length(unique(groups[,2]))){
    cells=row.names(subset(groups, groups[,1]==clust))
    subsetResults=resultsreadable[row.names(resultsreadable) %in% cells,]
    averagesReal[clust,]=apply(subsetResults,2,mean)
    rownames(averagesReal)[clust]=as.character(subset(groups, groups[,1]==clust)[1,2])
  }

  #create a table with average profiles of cell clusters and synthetic combinations
  allProfiles=rbind(averagesReal, synthProfiles)

  #this section determines the profile with the highest correlation to the given cell and determines if it is one of the doublet profiles
  for(cell in 1:nrow(isADoublet)){
    correlations=apply(allProfiles, 1, cor, resultsreadable[cell,])
    maxCorrelations=sort(correlations, decreasing = T)[1:2]
    maxCorrelation1=which(correlations==maxCorrelations[1])
    maxCorrelation2=which(correlations==maxCorrelations[2])
    #chosenCorrelation=min(maxCorrelation1, maxCorrelation2)
    chosenCorrelation=maxCorrelation1
    isADoublet[cell,1]=correlations[chosenCorrelation]
    correlatedCluster=row.names(allProfiles)[chosenCorrelation]
    isADoublet[cell,2]=correlatedCluster
    if(chosenCorrelation>length(unique(groups[,2]))){
    #if((chosenCorrelation>length(unique(groups[,2]))))&(maxCorrelation1-maxCorrelation2 >= 0.1)){
      isADoublet[cell,3]=TRUE
    }else{
      isADoublet[cell,3]=FALSE
    }
  }
  isADoublet[,4]=groups[,2]

  colnames(isADoublet)=c("Correlation","Cell_Types","isADoublet","Group_Cluster")

  return(list(isADoublet=isADoublet, resultsreadable=resultsreadable))

}
