#' Improved Seurat Pre Process
#'
#' This function improves the existing Seurat_Pre_Process function by working directly with a Seurat 3 or 4 object.
#' @param seuratObject Seurat object following a protocol such as https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
#' @param num_genes Number of genes for the top_n function. Default is 50.
#' @param write_files Save the output files to .txt format. Default is FALSE.
#' @param data_type Data slot to pull expression from in the Seurat object: counts or data. Default is "counts".
#' @return newExpressionFile - Seurat expression file in ICGS format (ICGS genes)
#' @return newFullExpressionFile - Seurat expression file in ICGS format (all genes)
#' @return newGroupsFile - Groups file ICGS format
#' @keywords Seurat
#' @export

Improved_Seurat_Pre_Process <- function(seuratObject, num_genes=50, write_files=FALSE, data_type="counts"){

  #Seurat version identification
  version = packageVersion("Seurat")

  #Upgrade Seurat object to the appropriate version for the version of Seurat you have installed
  seuratObject=UpdateSeuratObject(object = seuratObject)

  #extract expression file
  if(data_type=="counts"){
    expression=as.data.frame(seuratObject@assays[["RNA"]]@counts)
  }else if(data_type=="data"){
    expression=as.data.frame(seuratObject@assays[["RNA"]]@data)
  }else if(data_type=="scaled.data"){
    expression=as.data.frame(seuratObject@assays[["RNA"]]@scale.data)
  }

  #extract marker genes
  seuratObject.markers=FindAllMarkers(object = seuratObject, only.pos = TRUE, min.pct=0.25)
  if(version >= package_version(x = "3.9.9")){
    genes=seuratObject.markers %>% group_by(cluster) %>% top_n(n = num_genes, wt = avg_log2FC)
  }else if(version >= package_version(x = "3.0.0") && version < package_version(x = '3.9.9')){
    genes=seuratObject.markers %>% group_by(cluster) %>% top_n(n = num_genes, wt = avg_logFC)
  }else{
    print("This function only works with Seurat 3 or 4. Please update Seurat.")
  }

  #extract clusters
  clusters=as.data.frame(Idents(object = seuratObject))

  #Find and replace "-"
  colnames(expression)=gsub("-",".",colnames(expression))
  if(class(clusters[,1])=="character"){ #I believe this condition will never be reached but I am leaving it in as legacy code in case there is an edge case I'm missing!
    clusters[,1]=gsub("-",".",clusters[,1])
  }else{
    row.names(clusters)=gsub("-",".",row.names(clusters)) #added to fix bug
  }

  #Start cluster numbers at 1
  if(class(genes$cluster)=="factor"){
    if(min(as.numeric(as.character(genes$cluster)))==0){
      genes$cluster=as.numeric(as.character(genes$cluster))+1
      clusters[,1]=as.numeric(as.character(clusters[,1]))+1
    }else if(min(as.numeric(as.character(genes$cluster)))==1){
      genes$cluster=as.numeric(as.character(genes$cluster))
      clusters[,1]=as.numeric(as.character(clusters[,1]))
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



