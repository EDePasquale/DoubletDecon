#' Pseudo MarkerFinder
#'
#' This function uses t-tests to look for unique gene expression in each cluster.
#' @param groups Processed groups file from Clean_Up_Input.
#' @param data Processed data from Clean_Up_Input (or Remove_Cell_Cycle).
#' @param full_data2 Variable for passing the full gene expression matrix.
#' @return hallmarkTable - pseudo marker finder t statistics (gene by cluster).
#' @return hallmarkTable2 - pseudo marker finder gene assignments (which cluster has the highest t-stat for each gene).
#' @keywords Marker Finder ttest
#' @export

Pseudo_Marker_Finder<-function(groups, data, full_data2){

  #Save mean expression matrix  by cluster
  if(!is.null(full_data2)){
    centroids=data.frame(rep(NA,(nrow(full_data2)-1)))
    data2=apply(full_data2, 2,  as.numeric)
  }else{
    centroids=data.frame(rep(NA,(nrow(data)-1)))
    data2=apply(data, 2,  as.numeric)
  }

  clusters=length(unique(groups[,1]))

  for(cluster in 1:clusters){
    if(is.matrix(data2[2:nrow(data2),which(data2[1,]==cluster)])){
      centroids=cbind(centroids,apply(data2[2:nrow(data2),which(data2[1,]==cluster)],1,mean))
    }else{
      centroids=cbind(centroids,data2[2:nrow(data2),which(data2[1,]==cluster)])
    }
  }
  centroids=centroids[,-1]
  colnames(centroids)=unique(groups[,2])
  if(!is.null(full_data2)){
    row.names(centroids)=row.names(full_data2)[2:nrow(full_data2)]
  }else{
    row.names(centroids)=row.names(data)[2:nrow(data)]
  }

  if(!is.null(full_data2)){
    flag=TRUE
    sendingData=full_data2
  }else{
    flag=FALSE
    sendingData=data
  }

  new_table=as.data.frame(matrix(nrow=nrow(centroids),ncol=3))
  colnames(new_table)=c("Gene", "Cluster", "P-value")

  helper_PMF <-function(gene, centroids, sendingData, flag){
    data=sendingData

    temp=sort(centroids[gene, ], TRUE)[1:2]

    #Find highest and next highest cluster
    if(temp[1]==temp[2]){
      temp1=which(as.numeric(centroids[gene,])==as.numeric(temp[1]))[1]
      temp2=which(as.numeric(centroids[gene,])==as.numeric(temp[1]))[2]
    }else{
      temp1=min(which(as.numeric(centroids[gene,])==as.numeric(temp[1])))#use min here in case there are 2 that have the same value as it doen't matter which one I pick
      temp2=min(which(as.numeric(centroids[gene,])==as.numeric(temp[2])))
    }

    #Pull values each and do a t-test
    temp3=apply(data[,which(data[1,]==temp1 | data[1,]==temp2)],2,as.numeric)

    #Make sure that there is more than 1 cell in the cluster and save p-value for t-test
    if(ncol(as.data.frame(temp3[,which(temp3[1,]==temp1)])) >1 & ncol(as.data.frame(temp3[,which(temp3[1,]==temp2)])) >1){ #May still need this
      temp4=t.test(temp3[gene+1, which(temp3[1,]==temp1)], temp3[gene+1, which(temp3[1,]==temp2)])$p.value
    }else{
      temp4=NA
    }

    #Track progress
    setTxtProgressBar(pb, gene)
    # if(gene %% 100 == 0){
    #   print(paste0(gene, temp4))
    # }

    #If the p-value is significant, save the gene and cluster info
    if(!is.na(temp4) & temp4 <= 0.05){
      new_table[gene, 1]<<-row.names(centroids)[gene]
      new_table[gene, 2]<<-temp1
      new_table[gene, 3]<<-temp4
    }
  }

  #TODO: this is the slow part
  pb <- txtProgressBar(min = 1,
                        max = nrow(centroids),
                        initial = 0,
                        char = "=",
                        width = NA,
                        style = 3,
                        file = "")
  sapply(1:nrow(centroids), #For all genes
          helper_PMF, #Function
          centroids, #genes by centroids
          sendingData, #expression data (either full or ICGS/Seurat)
          flag) #TRUE if sending full data
  close(pb)

  new_table=new_table[complete.cases(new_table), ]

  return(new_table)
}
