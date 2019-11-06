#' Improved Seurat Pre Process
#'
#' This function improves the existing Seurat_Pre_Process function by working directly with a Seurat 3 object.
#' @param seuratObject Seurat object following a protocol such as https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
#' @param num_genes Number of genes for the top_n function. Default is 50.
#' @param write_files Save the output files to .txt format. Defauly is FALSE.
#' @return newExpressionFile - Seurat expression file in ICGS format (ICGS genes)
#' @return newFullExpressionFile - Seurat expression file in ICGS format (all genes)
#' @return newGroupsFile - Groups file ICGS format
#' @keywords Seurat
#' @export

Improved_Seurat_Pre_Process <- function(seuratObject, num_genes=50, write_files=FALSE){

  #load dplyr
  require(dplyr)
  
  #extract expression file
  expression=as.data.frame(seuratObject@assays[["RNA"]]@counts)
  
  #extract top 50 genes
  seuratObject.markers=FindAllMarkers(object = seuratObject, only.pos = TRUE, min.pct=0.25, logfc.threshold = 0.25)
  cluster1.markers=FindMarkers(object = seuratObject, ident.1 =0, thresh.use = 0.25, test.use = "roc", only.pos=TRUE)
  genes=seuratObject.markers %>% group_by(cluster) %>% top_n(n = num_genes, wt = avg_logFC)
  
  #extract clusters
  clusters=as.data.frame(Idents(object = seuratObject))

  #Find and replace "-"
  colnames(expression)=gsub("-",".",colnames(expression))
  clusters[,1]=gsub("-",".",clusters[,1])

  #Start cluster numbers at 1
  if(class(genes$cluster)=="factor"){
    if(min(as.numeric(as.character(genes$cluster)))==0){
      genes$cluster=as.numeric(as.character(genes$cluster))+1
      clusters[,1]=as.numeric(clusters[,1])+1
    }else if(min(as.numeric(as.character(genes$cluster)))==1){
      genes$cluster=as.numeric(as.character(genes$cluster))
      clusters[,1]=as.numeric(clusters[,1])
    }else{
      print("Unexpected cluster numbering scheme. Cluster numbers are expected to be continuous numbers starting from either 0 or 1. Please check conversion for correctness following this function.")
    }
  }
  

  #Reorder clusters file based on cluster number (this will be the new order for the cells)
  clusters2=clusters[order(clusters[,1]), , drop=FALSE]

  #Reorder columns
  expression=expression[row.names(clusters2)]

  #Reorder genes file based on gene cluster number (this will be the new order for the genes)
  genes2=genes[order(genes$cluster),]

  #Subset expression file for genes in the genes file
  allgenes=expression
  expression=expression[row.names(expression) %in% as.character(genes$gene),]

  #Reorder genes
  geneOrder=intersect(genes2$gene, as.character(row.names(expression)))
  expression=expression[match(geneOrder, row.names(expression)),]

  #Add columnn_clusters-flat
  allgenes=rbind(clusters2[,1], allgenes)
  expression=rbind(clusters2[,1], expression)
  row.names(allgenes)[1]="column_clusters-flat"
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
  
  if(write_files==TRUE){
    write.table(expression, "ICGS_expression.txt", sep="\t")
    write.table(allgenes, "ICGS_fullExpression.txt", sep="\t")
    write.table(groups, "ICGS_groups.txt", sep="\t", col.names = F)
  }

  return(list(newExpressionFile=expression,
              newFullExpressionFile=allgenes,
              newGroupsFile=groups))
}


