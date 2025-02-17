#' Penalty-Lift Analysis
#'
#' Penalty-Lift analysis for CATA variables, which is the difference between
#' the average hedonic response when CATA attribute is checked vs. the average
#' hedonic response when CATA attribute is not checked. 
#' @name plift
#' @usage plift(X, Y, digits = getOption("digits"), verbose = FALSE)
#' @param X  \emph{either} a matrix of CATA data with \eqn{I} consumers (rows)
#' and \eqn{J} products (columns) \emph{or} an array of CATA data with 
#' \eqn{I} consumers, \eqn{J} products, and \eqn{M} attributes.
#' @param Y  matrix of hedonic data with \eqn{I} consumers (rows)
#' and \code{J} products (columns)
#' @param digits for rounding
#' @param verbose set to \code{TRUE} to report counts and averages for checked
#' and not checked conditions (default: \code{FALSE})
#' @return Penalty lift per attribute, with counts and averages if \code{verbose} 
#' is \code{TRUE}.
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
#' plift(bread$cata[1:12,,1], bread$liking[1:12, ], digits = 3) 
#' 
#' # for the first 3  attributes with counts and averages
#' plift(bread$cata[1:12,,1:3], bread$liking[1:12, ], digits = 3, verbose = TRUE) 
plift <- function(X, Y, digits = getOption("digits"), verbose = FALSE){
  .pLiftAtt <- function(X, Y){
    if(length(X) != length(Y)){
      return(c(n1 = 0, n0 = 0, x1 = NA, x0 = NA, plift = NA))
    }
    if(length(X) == 0){
      return(c(n1 = 0, n0 = 0, x1 = NA, x0 = NA, plift = NA))
    }
    # now is the time to drop misssings
    xy <- matrix(c(X, Y), ncol=2)
    xy.drop.indx <- which(apply(xy, 1, function(.x){
      any(is.na(.x)) }))
    if(length(xy.drop.indx)>0){
      xy <- xy[-xy.drop.indx,]
    }
    if(nrow(xy)>0){
      indx1 <- which(xy[,1]==1)
      indx0 <- which(xy[,1]==0)
      x1 <- mean(xy[indx1, 2]) 
      x0 <- mean(xy[indx0, 2])
      return(c(n1 = length(indx1), n0 = length(indx0), x1 = x1, x0 = x0, plift = x1-x0))
    } else {
      return(c(n1 = 0, n0 = 0, x1 = NA, x0 = NA, plift = NA))
    }
  }
  if(!all.equal(dim(X)[1:2], dim(Y)[1:2])){
    # print("X and Y incompatible")
    # return(c(n1 = 0, n0 = 0, x1 = NA, x0 = NA, plift = NA))
  }
  if(length(dim(X)) == 2){
    X <- matrix(as.integer(X), nrow = nrow(X), ncol = ncol(X), dimnames=dimnames(X))
    if(!all(unique(c(X)) %in% 0:1)){
      # "X must be CATA data (coded 0/1)"
    }
    xdrop.indx <- which((rowSums(X) == ncol(X)) | (rowSums(X) == 0))
    # Check for flatliners will simply remove them rather than set them to missings
    # since X and Y are already 2-dimensional
    X[xdrop.indx,] <- NA 
    o <- .pLiftAtt(X, Y)
  } else if (length(dim(X)) == 3) {
    X <- array(as.integer(X), dim=dim(X), dimnames=dimnames(X))
    if(!all(unique(c(X)) %in% 0:1)){
      # "X must be CATA data (coded 0/1)"
    }
    # Set flatliners in X to NA. They will be removed later.
    Xsd <- apply(X, c(1,3), stats::sd)
    Xsd[Xsd == 0] <- NA
    Xsd[!is.na(Xsd)] <- 1
    Xna <- sweep(X, c(1,3), Xsd, "*")
    o <- apply(Xna, 3, function(x){ .pLiftAtt(x, Y) })
  }
  o <- t(as.matrix(o))
  o[,1:2] <- as.integer(round(o[,1:2], digits=0))
  o[,3:5] <- round(as.numeric(o[,3:5]), digits=digits)
  return(o)
}

