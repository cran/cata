#' Penalty-Lift Analysis
#'
#' Penalty-Lift analysis for CATA variables, which is the difference between
#' the average hedonic response when CATA attribute is checked vs. the average
#' hedonic response when CATA attribute is not checked. 
#' @name pLift
#' @aliases pLift
#' @usage pLift(X, Y)
#' @param X  \emph{either} a matrix of CATA data with \code{I} consumers (rows)
#' and \code{J} products (columns) \emph{or} an array of CATA data with 
#' \code{I} consumers, \code{J} products, and \code{M} attributes.
#' @param Y  matrix of hedonic data with \code{I} consumers (rows)
#' and \code{J} products (columns)
#' @return Penalty lift for the attribute if \code{X} is a matrix; otherwise,
#' penalty-lift for each attribute if \code{X} is a 3d array.
#' @export
#' @encoding UTF-8
#' @references 
#' Meyners, M., Castura, J.C., & Carr, B.T. (2013). Existing and new 
#' approaches for the analysis of CATA data. \emph{Food Quality and Preference}, 
#' 30, 309-319, \doi{10.1016/j.foodqual.2013.06.010}
#' @examples
#' data(bread)
#' 
#' # penalty lift, based only on the first 12 consumers
#' 
#' # for the first attribute ("Fresh")
#' pLift(bread$cata[1:12,,1], bread$liking[1:12, ]) 
#' 
#' # for the first 3  attributes
#' pLift(bread$cata[1:12,,1:3], bread$liking[1:12, ]) 
pLift <- function(X, Y){
  .pLiftAtt <- function(x, y){
    return(mean(c(y)[c(x)==1]) - mean(c(y)[c(x)==0]))
  }
  if(!all(unique(c(X)) %in% 0:1)){
    "X must be CATA data (coded 0/1)"
  }
  if(!all.equal(dim(X)[1:2], dim(Y)[1:2])){
    return("X and Y incompatible")
  }
  if(length(dim(X)) == 2){
    o <- .pLiftAtt(X, Y)
  } else if (length(dim(X)) == 3) {
    cy <- c(Y)
    o <- apply(X, 3, function(x){ .pLiftAtt(x, Y) })
  }
  return(o)
}

