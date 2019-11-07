#' ICGS2_to_ICGS1
#'
#' This function converts the new ICGS2 expression file format to the ICGS1 format. In making this change the cluster numbers (column 1 in the groups file) are retained and not the cluster names (column 2 in the groups file).
#' @param rawDataFile ICGS expression or counts file (in ICGS expression file format).
#' @param groupsFile ICGS groups file.
#' @param log_file_name used for saving run notes to log file
#' @return processed - data.frame of genes by samples with a row of cell clusters (column_clusters-flat) and a column of gene clusters (row_clusters-flat) when available.
#' @return groups - groups file with cell names matching the expression file.
#' @keywords ICGS
#' @export

ICGS2_to_ICGS1 <-function(rawDataFile, groupsFile, log_file_name=NULL){
  
  #pull ICGS2 header (or have it fed in)
  exp.header=read.table(rawDataFile, sep="\t", nrows=1, stringsAsFactors = F, row.names=1)
  
  #split away the useless groups names, keep only the cell names. save new column names for later use.
  exp.header=gsub(".*:", "", exp.header)
  
  #read in expression file (or have it fed in)
  exp.ICGS=read.table(rawDataFile, sep="\t", stringsAsFactors = F, row.names=1, header=T)
  
  #replace the column names with the saved ones
  colnames(exp.ICGS)=exp.header
  
  #read in the groups file (or have it fed in)
  groups.ICGS=read.table(groupsFile, sep="\t", stringsAsFactors = F, row.names=1, header=F)
  
  #replace the column-clusters_flat with the first column of groups
  exp.ICGS[1,2:ncol(exp.ICGS)]=groups.ICGS[,1]
  
  #split the row names into two vectors, the useless groups names and the gene names
  if(length(grep(":", rownames(exp.ICGS)[2]))==1){
    groups.names=gsub(":.*", "", row.names(exp.ICGS))
    groups.genes=gsub(".*:", "", row.names(exp.ICGS))
  }else{
    groups.names=exp.ICGS[,1]
    groups.genes=row.names(exp.ICGS)
  }

  #replace the row names with the gene names
  row.names(exp.ICGS)=groups.genes
  
  #match the useless groups names to the good groups names and put these groups into the row-clusters_flat column.
  exp.ICGS[2:nrow(exp.ICGS),1]=as.numeric(mapvalues(groups.names[2:length(groups.names)], unique(groups.names[2:length(groups.names)]), unique(groups.ICGS[,1])))
  
  #write to the log file
  if(!is.null(log_file_name)){
    cat("ICGS2 file formatted", file=log_file_name, append=TRUE, sep="\n")
    cat("ICGS2 file formatted", sep="\n")
  }
  
  #return the new expression and groups files
  return(list(rawData=exp.ICGS, groups=groups.ICGS))
}