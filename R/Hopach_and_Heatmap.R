#' Hopach and Heatmap
#'
#' This function performs hopach clustering and creates a heatmap and groups file.
#' @param partialData Subset of data, either doublets or nondoublets called by DeconRNASeq.
#' @param fullData Entire processed data from Clean_Up_Data.
#' @param filename Unique filename to be incorporated into the names of outputs from the functions.
#' @param groups Processed groups file from Clean_Up_Input.
#' @return partialGroups - new groups file for the partialData based on hopach clustering.
#' @keywords HOPACH heatmap
#' @export

Hopach_and_Heatmap<-function(partialData, fullData, filename, groups){

  gdist.cosangle<-distancematrix(partialData,d="cosangle")
  gho.cosangle<-hopach(partialData,dmat=gdist.cosangle,ord="own")
  makeoutput(partialData,gho.cosangle,file=filename)

  #from hopach
  p <- nrow(partialData)
  out <- row.names(partialData)
  out <- data.frame(UID = out)
  out <- data.frame(Index = (1:p), out)
  uclust <- sort(unique(gho.cosangle$clust$labels))
  clust <- NULL
  for (i in 1:length(uclust)) {
    clust[gho.cosangle$clust$label == uclust[i]] <- (i - 1)
  }
  out <- data.frame(out, Cluster.Number = clust)
  out <- data.frame(out, Cluster.Label = gho.cosangle$clustering$labels)
  out <- out[gho.cosangle$clustering$order, ]
  out <- data.frame(out, Cluster.Level.Order = (1:p))
  out <- out[order(out[, 1]), ]
  out <- data.frame(out, Final.Label = gho.cosangle$final$label)
  out <- out[gho.cosangle$final$order, ]
  out <- data.frame(out, Final.Level.Order = 1:p)
  #end from hopach
  test=fullData[,match(out$UID, colnames(fullData)) ]
  test=cbind(fullData[,1], test)
  colnames(test)[1]="row_clusters-flat"
  test[1,2:ncol(test)]=substring(out$Cluster.Label,1,1) #add new cluster labels
  #test[2:nrow(test),1]=

  #make heatmap of hopach clustering
  breaks=seq(0, #start point of color key
             10,  #end point of color key
             by=0.05) #length of sub-division
  mycol <- colorpanel(n=length(breaks)-1,low="black",high= "yellow") #heatmap colors
  if(is.na(test[2,1])){
    heatmap<-heatmap.2(data.matrix(test[2:nrow(test), 2:ncol(test)]), #the data matrix
                       Colv=FALSE, # No clustering of columns
                       Rowv = FALSE, #no clustering of rows
                       col=mycol, #colors used in heatmap
                       ColSideColors = as.color(Renumber(test[1,2:ncol(test)]), alpha=1, seed=4), #column color bar
                       #RowSideColors = rowColors, # row color bar
                       breaks=breaks, #color key details
                       trace="none", #no trace on map
                       na.rm=TRUE, #ignore missing values
                       margins = c(5,5), # size and layout of heatmap window
                       xlab = "Samples", #x axis title
                       ylab =  "Genes", # y axis title
                       main = filename) #main title
  }else{
    heatmap<-heatmap.2(data.matrix(test[2:nrow(test), 2:ncol(test)]), #the data matrix
                       Colv=FALSE, # No clustering of columns
                       Rowv = FALSE, #no clustering of rows
                       col=mycol, #colors used in heatmap
                       ColSideColors = as.color(Renumber(test[1,2:ncol(test)]), alpha=1, seed=4), #column color bar
                       RowSideColors = as.color(Renumber(test[2:nrow(test),1]), alpha=1, seed=2), # row color bar
                       breaks=breaks, #color key details
                       trace="none", #no trace on map
                       na.rm=TRUE, #ignore missing values
                       margins = c(5,5), # size and layout of heatmap window
                       xlab = "Samples", #x axis title
                       ylab =  "Genes", # y axis title
                       main = filename) #main title
  }


  partialGroups=cbind(Renumber(test[1,2:ncol(test)]),Renumber(test[1,2:ncol(test)]))
  row.names(partialGroups)=colnames(test[2:ncol(test)])
  colnames(partialGroups)=colnames(groups)

  return(partialGroups)

}
