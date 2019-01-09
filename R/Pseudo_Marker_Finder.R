#' Pseudo MarkerFinder
#'
#' This function uses ANOVA to look for unique gene expression in each possible doublet cluster.
#' @param groups Processed groups file from Clean_Up_Input.
#' @param data Processed data from Clean_Up_Input (or Remove_Cell_Cycle).
#' @param full_data2 cleaned full expression matrix from Clean_Up_Input.
#' @param downsample allows for downsampling of cells when using full expression matrix (use with large datasets), default is "none".
#' @param sample_num number of cells per cluster with downsampling with "even", percent of cluster with "prop".
#' @param min_uniq minimum number of unique genes required for a cluster to be rescued
#' @param log_file_name used for saving run notes to log file
#' @return new_table - non-doublet clusters, as determined by the "Remove" and "Rescue" steps.
#' @keywords Marker Finder ANOVA
#' @export

Pseudo_Marker_Finder<-function(groups, data, full_data2, downsample="none", sample_num=NULL, min_uniq=4, log_file_name=log_file_name){

  #Define possible doublet clusters
  doublets=cbind(unique(groups[grep("even|one|two", groups[,2]),1]),unique(groups[grep("even|one|two", groups[,2]),1]))

  #Deal with full data
  if(!is.null(full_data2)){
    data=full_data2
  }

  #Reduce gene list with full data file based on FC cutoffs and minimum expression values
  if(!is.null(full_data2)){
    matrix_data=apply(as.matrix(data[2:nrow(data), 2:ncol(data)]),2,as.numeric)
    keep.exprs1 <- rowSums(matrix_data)>0 #minimum expression values by gene, currently needs to be expressed at all in any cell
    f=as.factor(as.character(data[1,2:ncol(data)]))
    i <- split(1:ncol(matrix_data), f)
    x <- sapply(i, function(i){ #this function gives mean values for each gene per cluster
      if(!is.null(nrow(matrix_data[,i]))){
        rowMeans(matrix_data[,i])
      }else{
        as.numeric(matrix_data[,i])
      }
    })
    x=log2(x+1)
    j <- combn(levels(f), 2)
    ret<-x[,j[1,]]-x[,j[2,]] #vectorized log fold change operation
    colnames(ret)=paste(j[1,], j[2,], sep = '-')
    max_log_fc=apply(abs(ret),1,max)
    keep.exprs2 <- max_log_fc>=log2(1.5)
    data=data[c(1, (which(keep.exprs1==TRUE & keep.exprs2==TRUE)+1)),] #reduce to genes that pass both filters
  }

  #Deal with downsampling
  new_data=matrix(nrow=nrow(data), ncol=1)
  row.names(new_data)=row.names(data)
  new_data[,1]=data[,1]
  if(downsample=="even"){
    for(i in 1:length(unique(groups[,2]))){
      index_cells=which(data[1,]==i)
      if(length(index_cells)==1){
        cells=as.data.frame(data[,index_cells], stringsAsFactors=FALSE)
      }else{
        cells=data[,index_cells]
      }
      colnames(cells)=colnames(data)[index_cells]
      if(ncol(cells)>sample_num){
        new_cells=sample(cells, sample_num)
      }else{
        new_cells=cells
      }
      new_data=cbind(new_data,new_cells)
    }
  }else if(downsample=="prop"){
    for(i in 1:length(unique(groups[,2]))){
      index_cells=which(data[1,]==i)
      if(length(index_cells)==1){
        cells=as.data.frame(data[,index_cells], stringsAsFactors=FALSE)
      }else{
        cells=data[,index_cells]
      }
      colnames(cells)=colnames(data)[index_cells]
      new_cells=sample(cells, ceiling(ncol(cells)*sample_num))
      new_data=cbind(new_data,new_cells)
    }
  }else{
    new_data=data
  }

  #Non-doublet clusters vs doublet clusters
  doub_clusters=as.numeric(unique(doublets[,1]))
  non_doub_clusters=as.numeric(unique(groups[which(!(as.numeric(groups[,1]) %in% doub_clusters)),1]))

  #Create table to store p-values
  hallmarkTable=as.data.frame(matrix(ncol=length(doub_clusters), nrow=nrow(new_data)-1), row.names = row.names(new_data)[2:nrow(new_data)])

  #ANOVA testing
  helper_PMF_a <-function(ind, doub_clust_test, temp_data){
    gene_data=as.data.frame(cbind(as.numeric(temp_data[ind,]),as.numeric(temp_data[1,])))
    gene_data[,2]<-factor(gene_data[,2])
    colnames(gene_data)=c("Sabor", "Tipo")
    a1=aov(Sabor ~ Tipo, gene_data)
    if((!is.na(summary(a1)[[1]][["Pr(>F)"]][1])) & summary(a1)[[1]][["Pr(>F)"]][1]<0.05){ #if anova is significant, move to next step
      posthoc <- TukeyHSD(x=a1, conf.level=0.95) #do post hoc test
      doub_posthoc <- posthoc[["Tipo"]][,"p adj"][grep(doub_clust_test, names(posthoc[["Tipo"]][,"p adj"]))] #pull out all post hoc results containing the doublet cluster
      if(length(which(doub_posthoc<0.05))==(length(unique(gene_data$Tipo))-1)){ #if the doublet cluster is significantly different from all other clusters, move to next step
        mean_results=sort(tapply(gene_data$Sabor, gene_data$Tipo, mean), decreasing=TRUE) #get the mean for all clusters, and sort in decreasing order
        if(names(mean_results[1])==doub_clust_test){ #if the doublet cluster has the highest expression then call it as a non doublet cluster
          return(min(doub_posthoc))
        }
      }
    }
    return(NA)
  }

  #Helper function for the ANOVA testing
  helper_PMF_a_temp <- function(ind, doub_clust_test){
    gene_data=as.data.frame(cbind(as.numeric(temp_data[ind,]),as.numeric(temp_data[1,])))
    gene_data[,2]<-factor(gene_data[,2])
    colnames(gene_data)=c("Sabor", "Tipo")
    return(NA)
  }


  #Check for unique expression in each of the possible doublet clusters
  rows=2:nrow(new_data)
  new_data2=apply(as.matrix(new_data),2,as.numeric)
  clust <- makeCluster(detectCores()) #variable number of cores depending on the system
  for(newCluster in 1:length(doub_clusters)){
    print(paste0("Processing cluster ", newCluster, "/", length(doub_clusters), "..."))
    doub_clust_test=doub_clusters[newCluster] #which doublet cluster are we testing for unique gene expression
    temp_data=new_data2[,which(new_data2[1,] %in% non_doub_clusters | new_data2[1,]==doub_clust_test)] #pull all cells that are part of that doublet cluster or any non-doublet cluster
    hallmarkTable[,newCluster]=parSapply(clust, rows, helper_PMF_a, doub_clust_test, temp_data)
  }
  stopCluster(clust)

  #Determine which clusters are considered to have unique expression
  colnames(hallmarkTable)=as.character(doub_clusters)
  unique_genes_by_cluster=sapply(1:ncol(hallmarkTable), function(x) length(which(!is.na(hallmarkTable[,x]))))
  cat(paste0("Unique Genes By Cluster: ", unique_genes_by_cluster), file=log_file_name, append=TRUE, sep="\n")
  unique_rescued_clusters=colnames(hallmarkTable)[which(unique_genes_by_cluster>=min_uniq)] #min genes unique
  all_rescued_clusters=c(as.character(unique_rescued_clusters), as.character(non_doub_clusters))
  new_table=as.data.frame(matrix(ncol=2, nrow=length(all_rescued_clusters)))
  new_table[,2]=all_rescued_clusters

  return(new_table)

}
