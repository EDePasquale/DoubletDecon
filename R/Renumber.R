#' Renumber
#'
#' This function renumbers a vector to avoid non-consecutive cluster lables generated from cell cycle removal, hopach clustering, or non-standard raw data input.
#' @param x Vector to be renumbered
#' @return y - renumbered vector
#' @keywords renumber
#' @export

Renumber<-function(x){
  y=as.integer(as.factor(as.numeric(x)))
  return(y)
}
