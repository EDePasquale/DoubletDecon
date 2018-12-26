#' Clean Up Input
#'
#' This function cleans up the raw data and groups files. It replaces the column clusters in the expression file with the clusters in the groups file. It ensures that all expression files have a column_clusters-flat and row_clusters-flat (even when data not provided it adds NAs). It uses Renumber as necessary to avoid non-consecutive cluster identifications.
#' @param rawData ICGS expression or counts file (in ICGS expression file format).
#' @param groups ICGS groups file.
#' @param rowClusters Optional vector containing row cluster information (used by Recluster). Default is NULL.
#' @param log_file_name used for saving run notes to log file
#' @return processed - data.frame of genes by samples with a row of cell clusters (column_clusters-flat) and a column of gene clusters (row_clusters-flat) when available.
#' @return groups - groups file with cell names matching the expression file.
#' @keywords clean
#' @export

Clean_Up_Input<-function(rawData, groups, rowClusters=NULL, log_file_name){

  if(row.names(rawData)[1] %in% "column_clusters-flat" && (colnames(rawData)[1] %in% "row_clusters.flat" || colnames(rawData)[1] %in% "row_clusters-flat")){ #standard ICGS, contains column and row clusters

    rawDataStrip=rawData[-1,]
    processed=rbind(c(NA,groups[,1]),rawDataStrip) #replace with provided groups file
    row.names(processed)[1]="column_clusters-flat"

  }else if(row.names(rawData)[1] %in% "column_clusters-flat"){ #only has column clusters

    rawDataStrip=rawData[-1,]
    processed=rbind(groups[,1],rawDataStrip) #replace with provided groups file
    row.names(processed)[1]="column_clusters-flat"
    processed=cbind(rep(NA, nrow(processed)), processed)
    colnames(processed)[1]="row_clusters.flat"

  }else if(colnames(rawData)[1] %in% "row_clusters.flat" || colnames(rawData)[1] %in% "row_clusters-flat"){ #only has row clusters

    processed=rbind(c(NA,groups[,1]),rawData) #replace with provided groups file
    row.names(processed)[1]="column_clusters-flat"

  }else{ #has no column or row clusters

    processed=rbind(groups[,1],rawData) #replace with provided groups file
    row.names(processed)[1]="column_clusters-flat"
    processed=cbind(rep(NA, nrow(processed)), processed)
    colnames(processed)[1]="row_clusters.flat"

  }
  if(is.null(rowClusters)==FALSE){
    #attach the row clusters information
    processed[2:nrow(processed),1]=rowClusters
  }

  rownames(groups)=colnames(processed)[2:ncol(processed)] #because row and column names with special characters can cause problems in R

  #Renumber to avoid non-consecutive cluster identifications in the data and groups files
  groups[,1]=Renumber(groups[,1])
  processed[1,2:ncol(processed)]=Renumber(processed[1,2:ncol(processed)])
  if(is.na(processed[2,1])==FALSE){
    processed[2:nrow(processed),1]=Renumber(processed[2:nrow(processed),1])
  }

  cat(paste0(ncol(processed)-2, " samples after processing"), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0(nrow(processed)-2, " genes after processing"), file=log_file_name, append=TRUE, sep="\n")

  return(list(processed=processed, groups=groups))
}
