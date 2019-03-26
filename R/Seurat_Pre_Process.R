#' Seurat Pre Process
#'
#' This function uses t-tests to look for unique gene expression in each cluster.
#' @param expressionFile Normalized expression file from Seurat
#' @param genesFile Gene list from Seurat
#' @param clustersFile Cluster list from Seurat
#' @return newExpressionFile - Seurat expression file in ICGS format
#' @return newGroupsFile - Groups file ICGS format
#' @keywords Seurat
#' @export

Seurat_Pre_Process <- function(expressionFile, genesFile, clustersFile){

  #Read in data
  expression=read.table(expressionFile, sep="\t",header=T, row.names=1)
  genes=read.table(genesFile, sep="\t",header=T, row.names=1)
  clusters=read.table(clustersFile, sep="\t",header=T)

  #Find and replace "-"
  colnames(expression)=gsub("-",".",colnames(expression))
  clusters[,1]=gsub("-",".",clusters[,1])

  #Start cluster numbers at 1
  genes$cluster=genes$cluster+1
  clusters$x=clusters$x+1

  #Reorder clusters file based on cluster number (this will be the new order for the cells)
  clusters2=clusters[order(clusters$x),]

  #Reorder columns
  expression=expression[as.character(clusters2$X)]

  #Reorder genes file based on gene cluster number (this will be the new order for the genes)
  genes2=genes[order(genes$cluster),]

  #Subset expression file for genes in the genes file
  allgenes=expression
  expression=expression[row.names(expression) %in% as.character(genes$gene),]

  #Reorder genes
  geneOrder=intersect(genes2$gene, as.character(row.names(expression)))
  expression=expression[match(geneOrder, row.names(expression)),]

  #Add columnn_clusters-flat
  allgenes=rbind(clusters2[,2], allgenes)
  expression=rbind(clusters2[,2], expression)
  row.names(genes2)[1]="column_clusters-flat"
  row.names(expression)[1]="column_clusters-flat"

  #Add row_clusters-flat
  genes3=genes2[match(geneOrder, genes2$gene),]
  rowToAdd=c(NA, genes3$cluster)
  rowToAdd2=rep(NA, nrow(allgenes))
  expression=cbind(rowToAdd, expression)
  allgenes=cbind(rowToAdd2, allgenes)
  colnames(expression)[1]="row_clusters-flat"
  colnames(allgenes)[1]="row_clusters-flat"

  #Make groups file
  groups=cbind(as.numeric(expression[1,2:ncol(expression)]), as.numeric(expression[1,2:ncol(expression)]))
  row.names(groups)=as.character(colnames(expression)[2:ncol(expression)])

  return(list(newExpressionFile=expression,
              newFullExpressionFile=allgenes,
              newGroupsFile=groups))
}


