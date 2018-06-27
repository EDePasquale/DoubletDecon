#' Seurat Pre Process
#'
#' This function uses t-tests to look for unique gene expression in each cluster.
#' @param expressionFile
#' @param genesFile
#' @param clustersFile
#' @return newExpressionFile - Seurat expression file in ICGS format
#' @return newGroupsFile - Groups file ICGS format
#' @keywords Seurat
#' @export

Seurat_Pre_Process <- function(expressionFile, genesFile, clustersFile){

  #Read in data
  expression=read.table(expressionFile, sep="\t",header=T, row.names=1)
  genes=read.table(genesFile, sep="\t",header=T, row.names=1)
  clusters=read.table(clustersFile, sep="\t",header=T)

  #Start cluster numbers at 1
  genes$cluster=genes$cluster+1
  clusters$x=clusters$x+1

  #Reorder clusters file based on cluster number (this will be the new order for the cells)
  clusters2=clusters[order(clusters$x),]

  #Reorder columns
  expression2=expression[clusters2$X]

  #Reorder genes file based on gene cluster number (this will be the new order for the genes)
  genes2=genes[order(genes$cluster),]

  #Subset expression file for genes in the genes file
  expression3=expression2[row.names(expression2) %in% as.character(genes$gene),]

  #Reorder genes
  geneOrder=intersect(genes2$gene, as.character(row.names(expression3)))
  expression4=expression3[match(geneOrder, row.names(expression3)),]

  #Add columnn_clusters-flat
  expression5=rbind(clusters2[,2], expression4)
  row.names(expression5)[1]="column_clusters-flat"

  #Add row_clusters-flat
  genes3=genes2[match(geneOrder, genes2$gene),]
  rowToAdd=c(NA, genes3$cluster)
  expression6=cbind(rowToAdd, expression5)
  colnames(expression6)[1]="row_clusters-flat"

  #Make groups file
  groups=cbind(as.numeric(expression6[1,2:ncol(expression6)]), as.numeric(expression6[1,2:ncol(expression6)]))
  row.names(groups)=as.character(colnames(expression6)[2:ncol(expression6)])

  return(list(newExpressionFile=expression6,
              newGroupsFile=groups))
}


