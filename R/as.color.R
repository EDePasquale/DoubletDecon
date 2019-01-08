#' Create a unique set of colors that matched the unique set of inputs from a single vector.
#' Original Author: silas.tittes - at- gmail.com
#' Removed from CRAN
#'
#' Provided a vector of any class, the function  will return a set of colors
#' with the same length and order as the input vector. A random seed is used
#' and will be printed in order to save desirable color palettes.
#' Color transparency can also be modified.
#'
#' @param x A vector of any class and length. \code{x}
#' @param alpha A numeric value between 0 and 1 to control color transparency. Default is value is 0.5. \code{alpha}
#' @param seed A value to seed the random colors. If not provided, a new palette of random colors will be generated with each call to the function.
#'
#' @return output a vector of hexadecimal colors of length x will be returned.
#'
#' @export
#'
#' @examples
#' #simple data frame with factors
#' n <- 100
#' f <- 5
#' x <- sort(rnorm(n, mean = 0, sd = 50)) + rnorm(n, mean = 0, sd = 30)
#' fact <- rep(letters[1:5], each=n/f)
#' #call to as.color, with char vector
#' colz <- as.color(fact, alpha = 1)
#' plotx <- as.integer(as.factor(fact))
#' plot( jitter(plotx), x, col=colz, pch=19)


as.color <- function(x, alpha, seed){

  #set random seed
  if(missing(seed)){
    lSeed <- round(runif(n = 1, min = 1, max = .Machine$integer.max), 0)
    print(paste("colors chosen using", lSeed, "as random seed", sep = " "))
    set.seed(lSeed)
  } else {
    set.seed(seed)
  }

  #set color alpha, or default to 0.5
  if(missing(alpha)){
    alpha <- 0.5
  }

  #get the number of colors needed
  numcols <- unique(x)
  colseq <- seq(0, 1, length.out=255)

  #check for number of input colors
  if( numcols > 255^3){
    warning("The input vector has more unique values than possible colors,
	    some colors will be recycled.")
  }

  #create colors
  r <- sample(colseq, length(numcols), replace = F)
  g <- sample(colseq, length(numcols), replace = F)
  b <- sample(colseq, length(numcols), replace = F)
  colz <- rgb(red = r, green = g, blue = b, alpha = alpha)

  #assign colors to input values
  colvec <- rep(NA, length(x))
  for(i in 1:length(numcols)){
    colvec[x %in% numcols[i] ] <- colz[i]
  }
  return(colvec)
}
