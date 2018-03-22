#' Pseudo MarkerFinder
#'
#' This function uses t-tests to look for unique gene expression in each cluster.
#' @param groups Processed groups file from Clean_Up_Input.
#' @param data Processed data from Clean_Up_Input (or Remove_Cell_Cycle).
#' @return hallmarkTable - pseudo marker finder t statistics (gene by cluster).
#' @return hallmarkTable2 - pseudo marker finder gene assignments (which cluster has the highest t-stat for each gene).
#' @keywords Marker Finder ttest
#' @export

Pseudo_Marker_Finder<-function(groups, data){

  #Create table to store
  hallmarkTable=as.data.frame(matrix(ncol=length(unique(groups[,1])), nrow=nrow(data)-1), row.names = row.names(data)[2:nrow(data)])

  #T-test for each gene
  x=data.matrix(data[2:nrow(data),2:ncol(data)])
  for(newCluster in 1:length(unique(groups[,1]))){
    fac=as.numeric(colnames(data)[2:ncol(data)] %in% row.names(subset(groups, groups[,1]==newCluster)))
    temp=mt.teststat(x,fac,test="t")
    # pval=2*(1-pt(abs(temp),df=(ncol(data)-1)))
    # pval=2*(1-pt(abs(temp),df=73.442))
    # #print(pval[1])
    # a=x[,fac==1]
    # b=x[,fac==0]
    # z=as.numeric(a[newCluster,2:ncol(a)])
    # y=as.numeric(b[newCluster,2:ncol(b)])
    # temp2=t.test(z,y)
    # print(temp[newCluster])
    # #print(temp2$statistic)
    # temp4=cbind(x[1,], fac)
    # colnames(temp4)=c("extra", "group")
    # temp4=as.data.frame(temp4)
    # temp5=with(temp4, t.test(extra[group == 1], extra[group == 0]))
    # print(temp5)
    if(is.nan(temp[1])){
      hallmarkTable[,newCluster]=rep(0,length(temp))
    }else{
      hallmarkTable[,newCluster]=temp
    }
  }

  #For each gene, which cluster has the strongest "claim"
  hallmarkTable2=as.data.frame(matrix(ncol=1, nrow=nrow(hallmarkTable)), row.names=row.names(hallmarkTable))
  for(gene in 1:nrow(hallmarkTable)){
    hallmarkTable2[gene,1]=which(hallmarkTable[gene,]==max(hallmarkTable[gene,]))
  }

  return(list(hallmarkTable=hallmarkTable, hallmarkTable2=hallmarkTable2))

}
