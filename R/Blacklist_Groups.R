#' Blacklist Groups
#'
#' This function is used for the calculation of blacklisted medoids and groups. It calculates medoids. It performs medoid correlations. It creates the binary correlation table between medoids, otherwise known as the blacklist. It makes a blacklist heatmap. It makes a blacklist heatmap. It makes new medoids based on blacklisted combined clusters. It makes a new groups file based on the blacklisted combined clusters.
#' @param data  Processed data from Clean_Up_Input (or Remove_Cell_Cycle)
#' @param groups Processed groups file from Clean_Up_Input
#' @param rhop x in mean+x*SD to determine upper cutoff for correlation in the blacklist. Default is 1.
#' @param centroid_flag use centroids if data is sparse
#' @param log_file_name used for saving run notes to log file
#' @return newMedoids - new medoids data.frame for the new combined blacklisted clusters.
#' @return newGroups - new groups file containing cluster assignment based on new combined blacklisted clusters.
#' @keywords blacklist correlation
#' @export

Blacklist_Groups<-function(data, groups, rhop, centroid_flag, log_file_name){

  #Step 1: calculate medoids
  medoids=data.frame(rep(NA,nrow(data)-1))
  clusters=length(unique(groups[,1]))
  for(cluster in 1:clusters){
    if(centroid_flag==TRUE){
      medoids=cbind(medoids,apply(data[2:nrow(data),which(data[1,]==cluster)],1,mean))
    }else{
      medoids=cbind(medoids,apply(data[2:nrow(data),which(data[1,]==cluster)],1,median))
    }
  }
  medoids=medoids[,-1]
  colnames(medoids)=unique(groups[,2])

  #Step 2: medoid correlations
  cormedoids=cor(medoids, method="pearson")

  #Step 3: create blacklist (binary correlation table between medoids)
  blacklist=data.frame(matrix(ncol=ncol(medoids),nrow=ncol(medoids)))
  cutoff=mean(cormedoids)+(rhop)*sd(cormedoids) #based on mean + 1SD * user provided multiplier
  for(rrow in 1:nrow(cormedoids)){
    for(ccol in 1:ncol(cormedoids)){
      if(cormedoids[rrow,ccol]>cutoff){
        blacklist[rrow,ccol]=1
      }else{
        blacklist[rrow,ccol]=0
      }
    }
  }
  blacklist_original_order=blacklist #keep original order with no names so I can pull the correct clusters out when trying to make medoids
  row.names(blacklist)=row.names(cormedoids)
  colnames(blacklist)=row.names(cormedoids)

  #Step 4: make blacklist heatmap
  BLheatmap=heatmap.2(as.matrix(blacklist_original_order),
                      Colv=TRUE, # clustering of columns
                      Rowv=TRUE, # clustering of rows
                      xlab = "Cell Types", #x axis title
                      ylab =  "Cell Types", #y axis title
                      trace="none",
                      main = "Cluster Merge") #main title
  blacklist_original_order=blacklist_original_order[BLheatmap$rowInd,BLheatmap$colInd]
  blacklist=blacklist[BLheatmap$rowInd,BLheatmap$colInd]

  #Step 5: make new medoids
  blacklistCluster=try(mcl(blacklist, addLoops=FALSE)$Cluster)
  if(class(blacklistCluster) == "try-error"){
    print("Unable to perform mcl function for blacklist clustering, please try a different rhop.")
    stop()
  }
  i=-1 #if the cluster is 0 (meaning no combining) change the name of the cluster to make it unique
  for(cluster in 1:length(blacklistCluster)){
    if(blacklistCluster[cluster]==0){
      blacklistCluster[cluster]=i
      i=i-1
    }
  }
  uniquelist=unique(blacklistCluster) #list of unique clusters
  nunique=length(uniquelist) #number of unique clusters
  newMedoids=data.frame(matrix(ncol=nunique, nrow=(nrow(data)-1)))
  for(cluster in 1:nunique){ #for each new cluster, assign medoid to new data.frame
    temp=which(blacklistCluster==uniquelist[cluster])
    if(centroid_flag==TRUE){
      newMedoids[,cluster]=apply(data[2:nrow(data),(data[1,] %in% rownames(blacklist_original_order)[temp])],1,mean) #this is where the original order is critical
    }else{
      newMedoids[,cluster]=apply(data[2:nrow(data),(data[1,] %in% rownames(blacklist_original_order)[temp])],1,median) #this is where the original order is critical
    }
    colnames(newMedoids)[cluster]=paste(rownames(blacklist)[temp], collapse="-") #need to give the columns meaningful names (combination names of combined clusters)

  }
  row.names(newMedoids)=row.names(data)[2:nrow(data)]

  #Step 6: make new combined groups file
  newGroups=data.frame(matrix(ncol=ncol(groups)+1, nrow=1))
  for(cluster in 1:nunique){
    temp=which(blacklistCluster==uniquelist[cluster])
    temp1.5=as.integer(row.names(blacklist_original_order)[temp])
    temp2=groups[groups[,1] %in% temp1.5,]
    temp3=cbind(temp2,rep(paste(rownames(blacklist)[temp], collapse="-"), nrow(temp2)))
    colnames(temp3)=colnames(newGroups)
    newGroups=rbind(newGroups, temp3)
  }
  newGroups=newGroups[-1,]
  newGroups[,2]=as.integer(as.factor(newGroups[,3]))
  newGroups=newGroups[,-1]

  cat(paste0("New blacklisted clusters: ", paste(unique(newGroups[,2]), sep="' '", collapse=", ")), file=log_file_name, append=TRUE, sep="\n")


  return(list(newMedoids=newMedoids, newGroups=newGroups))
}
