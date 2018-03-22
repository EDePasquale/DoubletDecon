#' Is A Doublet 2
#'
#' This function uses deconvolution analysis (DeconRNASeq) to evaluate each cell for equal contribution from blacklisted clusters.
#' @param data Processed data from Clean_Up_Input (or Remove_Cell_Cycle).
#' @param newMedoids New combined medoids from Blacklist_groups.
#' @param groups Processed groups file from Clean_Up_Input.
#' @param synthProfiles Average profiles of synthetic doublets from Synthetic_Doublets.
#' @return isADoublet - data.frame with each cell as a row and whether it is called a doublet by deconvolution analysis.
#' @return resultsreadable - data.frame with results of deconvolution analysis (cell by cluster) in percentages.
#' @keywords doublet deconvolution decon
#' @export

Is_A_Doublet2<-function(data, newMedoids, groups, synthProfiles){

  doubletCutoff=30

  #create data frame to store doublets table
  isADoublet=data.frame(matrix(ncol=6,nrow=(ncol(data)-1)))
  rownames(isADoublet)=colnames(data)[2:ncol(data)]
  rownames(newMedoids)=rownames(data)[2:nrow(data)]

  #run DeconRNASeq with new medoids and data
  results=DeconRNASeq(data[2:nrow(data), 2:ncol(data)], newMedoids)
  resultsreadable=round(results$out.all*100,2)
  rownames(resultsreadable)=rownames(isADoublet) #make an easily readable results table

  #this section determines if the biggest 2 contributors are larger than the doublet cutoff (ex: 30%)
  for(cell in 1:nrow(isADoublet)){
    meansSort=sort(resultsreadable[cell,], decreasing=T)
    if(meansSort[1]>doubletCutoff && meansSort[2]>doubletCutoff){
      if(meansSort[1]==meansSort[2]){ #to catch the rare occurance that two medoids equally contribute to a cell
        isADoublet[cell,1]=TRUE
        temp=which(resultsreadable[cell,]==meansSort[[1]])
        isADoublet[cell,2]=colnames(resultsreadable)[temp[1]]
        isADoublet[cell,3]=meansSort[[1]]
        isADoublet[cell,4]=colnames(resultsreadable)[temp[2]]
        isADoublet[cell,5]=meansSort[[2]]
        if(isADoublet[cell,2]<isADoublet[cell,4]){
          isADoublet[cell,6]=paste0(isADoublet[cell,2], "_", isADoublet[cell,4])
        }else{
          isADoublet[cell,6]=paste0(isADoublet[cell,4], "_", isADoublet[cell,2])
        }
      }else{
        isADoublet[cell,1]=TRUE
        isADoublet[cell,2]=colnames(resultsreadable)[which(resultsreadable[cell,]==meansSort[[1]])]
        isADoublet[cell,3]=meansSort[[1]]
        isADoublet[cell,4]=colnames(resultsreadable)[which(resultsreadable[cell,]==meansSort[[2]])]
        isADoublet[cell,5]=meansSort[[2]]
        if(isADoublet[cell,2]<isADoublet[cell,4]){
          isADoublet[cell,6]=paste0(isADoublet[cell,2], "_", isADoublet[cell,4])
        }else{
          isADoublet[cell,6]=paste0(isADoublet[cell,4], "_", isADoublet[cell,2])
        }
      }
    }else{
      isADoublet[cell,1]=FALSE
      isADoublet[cell,2]=NA
      isADoublet[cell,3]=NA
      isADoublet[cell,4]=NA
      isADoublet[cell,5]=NA
      isADoublet[cell,6]=NA
    }
  }
  isADoublet=cbind(isADoublet, groups[,2])
  colnames(isADoublet)=c("isADoublet","cellType1","percentage1","cellType2","percentage2", "cluster", "Cell_Types")

  return(list(isADoublet=isADoublet, resultsreadable=resultsreadable))

}
