#' Calculate within-cluster homogeneity
#'
#' Within a group of \code{N} consumers, the Homogeneity index lies between
#' \code{1/N} (no homogeneity) to \code{1} (perfect homogeneity). 
#' @name homogeneity
#' @aliases homogeneity
#' @usage homogeneity(X)
#' @param X three-way array; the \code{I, J, M} array has \code{I}
#' assessors, \code{J} products, code{M} attributes where CATA data have values 
#' \code{0} (not checked) and \code{1} (checked)
#' @return homogeneity index
#' @export
#' @encoding UTF-8
#' @references Llobell, F., Cariou, V., Vigneau, E., Labenne, A., & Qannari, 
#' E. M. (2019). A new approach for the analysis of data and the clustering of  
#' subjects in a CATA experiment. \emph{Food Quality and Preference}, 72, 31-39, 
#' \doi{10.1016/j.foodqual.2018.09.006}
#' 
#' @examples
#' data(bread)
#' 
#' # homogeneity index for the first 7 consumers on the first 6 attributes
#' homogeneity(bread$cata[1:7,,1:6])
homogeneity <- function(X){
  nI <- dim(X)[1]
  nJ <- dim(X)[2]
  nM <- dim(X)[3]
  i.norm <- array(0, dim = c(nJ, nM, nI))
  if(any(apply(X, 1, sum)==0)) return("Some assessors give no responses.")
  for (i in 1:nI) {
    i.cata <- as.matrix(X[i,,])
    i.norm[, , i] <- i.cata/sqrt(sum(i.cata == 1))
  }
  S <- matrix(0, nI, nI)
  diag(S) <- rep(1, nI)
  for (i1 in 1:(nI - 1)) {
    for (i2 in (i1 + 1):nI) {
      S[i2, i1] <- S[i1, i2] <- sum(diag(tcrossprod(i.norm[, , i1], 
                                                    i.norm[, , i2])))
    }
  }
  this.lambda <- 1 
  if(nI > 1){
    this.lambda  <-  svd(S)$d[1]
  }
  return(100*this.lambda/nI)
}

#' Evaluate Quality of Cluster Analysis Solution
#'
#' Evaluate the quality of cluster analysis solutions using measures related to
#' within-cluster product discrimination, between-cluster non-redundancy,
#' overall diversity (coverage), average RV, sensory differentiation retained,
#' and within-cluster homogeneity. 
#' @name evaluateClusterQuality
#' @aliases evaluateClusterQuality
#' @usage evaluateClusterQuality(X, M, alpha = .05, M.order = NULL, 
#' quiet = FALSE, digits = getOption("digits"), ...)
#' @param X three-way array; the \code{I, J, M} array has \code{I}
#' assessors, \code{J} products, code{M} attributes where CATA data have values 
#' \code{0} (not checked) and \code{1} (checked)
#' @param M  cluster memberships
#' @param alpha significance level to be used for two-tailed tests
#' @param M.order can be used to change the cluster numbers (e.g. to label 
#' cluster 1 as cluster 2 and vice versa); defaults to \code{NULL}
#' @param quiet if \code{FALSE} (default) then it prints information quality
#' measures; if \code{TRUE} then returns results without printing
#' @param digits significant digits (to display)
#' @param ... other parameters for \code{\link[base]{print.default}} (if 
#' \code{quiet = TRUE}).
#' @return A list containing cluster analysis quality measures: 
#' \itemize{
#'  \item{\code{$solution} : 
#'    \itemize{
#'    \item{\code{Pct.b} = percentage of the total sensory differentiation 
#'    retained in the solution}
#'    \item{\code{min(NR)} = smallest observed between-cluster non-redundancy}
#'    \item{\code{Div_G} = overall diversity (coverage)}
#'    \item{\code{H_G} = overall homogeneity (weighted average of within-cluster
#'    homogeneity indices)}
#'    \item{\code{avRV} = average RV coefficient for all between-cluster 
#'    comparisons}}}
#'  \item{\code{$clusters} : 
#'    \itemize{
#'    \item{\code{ng} = number of cluster members}
#'    \item{\code{bg} = sensory differentiation retained in cluster}
#'    \item{\code{xbarg} = average citation rate in cluster}
#'    \item{\code{Hg} = homogeneity index within cluster (see 
#'    \code{\link[cata]{homogeneity}})}
#'    \item{\code{Dg} = within-cluster product discrimination}}}
#'  \item{\code{$nonredundancy.clusterpairs} : 
#'    \itemize{
#'    \item{square data frame showing non-redundancy for each pair of clusters
#'    (low values indicate high redundancy)}}}
#'  \item{\code{$rv.clusterpairs} : 
#'    \itemize{
#'    \item{square data frame with RV coefficient for each pair of clusters
#'    (high values indicate higher similarity in product configurations)}}}}
#' @export
#' @encoding UTF-8
#' @seealso \code{\link[cata]{homogeneity}}
#' @references Castura, J.C., Meyners, M., Varela, P., & Næs, T. (2022). 
#' Clustering consumers based on product discrimination in check-all-that-apply 
#' (CATA) data. \emph{Food Quality and Preference}, 104564. 
#' \doi{10.1016/j.foodqual.2022.104564}.
#' 
#' @examples
#' data(bread)
#' evaluateClusterQuality(bread$cata[1:14,,1:6], M = rep(1:2, each = 7))
evaluateClusterQuality <- function(X, M, alpha = .05, M.order = NULL, 
                       quiet = FALSE, digits = getOption("digits"), ...){
  # start of functions
  # calculate Within-cluster product discrimination (D_{g})
  .withinGroupDiscrim.signif <- function(b, c, 
                                         alternative = "two.sided", alpha=.05){
    if(b+c == 0) return(NA)
    x <- b
    if(alternative=="two.sided"){
      x <- min(b,c)
    } 
    return(stats::binom.test(x, b+c, alternative = alternative)$p.value < alpha)
  }
  # end of functions
  
  X.bc <- barray(X)
  gID <- sort(unique(M))
  if(!is.null(M.order) & length(M.order)==length(gID)) gID <- gID[M.order]
  G <- length(gID)
  
  # *** Results per cluster *** 
  gMat <- matrix(0, nrow=G, ncol=5, 
                 dimnames = list(1:G, c("ng", "bg", "xbarg", "Hg", "Dg")))
  # matrix for two-tailed p-values
  # probably this can be refactored to compare p-values from mcnemarQ() to alpha
  g.pp.M.two1tail <- array(0, c(G, dim(X.bc)[c(2,4)], 2),
                           dimnames = list(as.character(gID), 
                                           NULL, 
                                           dimnames(X.bc)[[4]],
                                           c("lower.tail", "upper.tail")))
  for(g in 1:G){
    Mg <- which(M == gID[g])
    g.pp.M.two1tail[g,,,1] <- mapply(.withinGroupDiscrim.signif,
                                     b = c(apply(X.bc[Mg, , 1, ], 2:3, sum)), 
                                     c = c(apply(X.bc[Mg, , 2, ], 2:3, sum)),
                                     MoreArgs = list(alternative = "l", 
                                                     alpha = alpha/2))
    g.pp.M.two1tail[g,,,2] <- mapply(.withinGroupDiscrim.signif,
                                     b = c(apply(X.bc[Mg, , 1, ], 2:3, sum)), 
                                     c = c(apply(X.bc[Mg, , 2, ], 2:3, sum)),
                                     MoreArgs = list(alternative = "g", 
                                                     alpha = alpha/2))
    gMat[g,1] <- length(Mg)
    gMat[g,2] <- getb(X.bc[Mg,,1,], X.bc[Mg,,2,])
    gMat[g,3] <- mean(X[Mg,,])
    # get homogeneity
    gMat[g,4] <- homogeneity(X[Mg,,])
    gMat[g,5] <- sum(g.pp.M.two1tail[g,,,], na.rm=TRUE)/
      prod(dim(g.pp.M.two1tail[1,,,1]))
  }
  
  # *** Non-redundancy results per pair of clusters *** 
  ggMat <- matrix(NA, nrow=G, ncol=G)
  
  # *** RV & Salton cosine results per pair of clusters *** 
  ggRVMat <- matrix(NA, nrow=G, ncol=G)
  #ggSaltonMat <- matrix(NA, nrow=G, ncol=G)
  if(G>1){
    for(gg1 in 2:nrow(ggMat)){
      g1.data <- toMatrix(X[which(M==gID[gg1]),,])
      g1.aggr <- stats::aggregate(g1.data[,-c(1,2)], 
                           list(g1.data[,2]), sum)[,-1]
      g1.prop <- as.matrix(g1.aggr / length(which(M==gID[gg1])))
      
      for(gg2 in 1:(gg1-1)){
        gg.diff.lower <- g.pp.M.two1tail[gg1,,,1] - g.pp.M.two1tail[gg2,,,1]
        gg.diff.upper <- g.pp.M.two1tail[gg1,,,2] - g.pp.M.two1tail[gg2,,,2]
        # above we are just getting rid of cases where both are significant as
        # it implies "no group-to-group difference"
        # gg.diff.lower and gg.diff.upper can have values in -1, 0, or 1
        # and their sum can be in -2, -1, 0, 1, 2
        # it is their absolute value that is of interest
        # Denominator for NR is 2*J*(J-1)*M 
        this.num <- sum(abs(gg.diff.lower), abs(gg.diff.upper), na.rm = TRUE)
        this.denom <- prod(dim(g.pp.M.two1tail)[-1])
        ggMat[gg1, gg2] <- 100*this.num/this.denom
        # RV
        g2.data <- toMatrix(X[which(M==gID[gg2]),,])
        g2.aggr <- stats::aggregate(g2.data[,-c(1,2)], 
                             list(g2.data[,2]), sum)[,-1]
        g2.prop <- as.matrix(g2.aggr / length(which(M==gID[gg2])))
        ggRVMat[gg1, gg2] <- rv.coef(g1.prop, g2.prop)
        # Salton
        #ggSaltonMat[gg1, gg2] <- salton(g1.prop, g2.prop)
      }
    }
  }
  # *** Solution results *** 
  Mat <- matrix(NA, nrow=1, ncol=5, dimnames = list(
    "Solution", c("Pct.b", "min(NR)", "Div_G", "H_G", "avRV")))
  Mat[1,1] <- 100*sum(gMat[,2])/sum(X.bc)
  Mat[1,2] <- ifelse(G>1, min(ggMat, na.rm=TRUE), NA)
  Mat[1,3] <- 100*sum(apply(g.pp.M.two1tail, c(2,3,4), 
                            function(x) (sum(x, na.rm = TRUE)>0)*1))/
    prod(dim(g.pp.M.two1tail)[-1])
  Mat[1,4] <- sum(gMat[,1] * 0.01*gMat[,4])/sum(gMat[,1])*100
  if(G>1){
    Mat[1,5] <- mean(ggRVMat[lower.tri((ggRVMat))], na.rm=TRUE)
  }
  
  # apply rounding to pairwise tables
  reportNR <- as.data.frame(ggMat)
  diag(reportNR) <- 0
  reportRV <- as.data.frame(ggRVMat)
  diag(reportRV) <- 1
  # reportSalton <- as.data.frame(ggSaltonMat)
  # diag(reportSalton) <- 1
  colnames(reportNR) <- colnames(reportRV) <- 
    rownames(reportNR) <- rownames(reportRV) <- 
    #rownames(reportSalton) <- rownames(reportSalton) <- 
    paste0("g",1:G)
  
  res <- list(
    solution = Mat, 
    clusters = gMat, 
    nonredundancy.clusterpairs = reportNR,
    rv.clusterpairs = reportRV #,
    #salton.clusterpairs = reportSalton
    )
  class(res) <- "cataQuality"
  if(!quiet){
    .unquote <- function(x, ...){ 
      return(print(x, ..., quote = FALSE))
    }
    o <- res
    # cleanup the look of the pairwise tables for printing
    o$solution <- signif(res$solution, digits = digits)
    o$clusters <- signif(res$clusters, digits = digits)
    o$nonredundancy.clusterpairs <- 
      signif(res$nonredundancy.clusterpairs, digits = digits)
    o$rv.clusterpairs <- 
      signif(res$rv.clusterpairs, digits = digits)
    o$nonredundancy.clusterpairs[
      upper.tri(o$nonredundancy.clusterpairs)] <- ""
    o$rv.clusterpairs[
      upper.tri(o$rv.clusterpairs)] <- ""
    # o$salton.clusterpairs[
    #   upper.tri(o$salton.clusterpairs)] <- ""
    cat(paste(c("",
                "Results", 
                "-------", ""), 
              collapse = '\n'))
    print(as.data.frame(o$solution, row.names = ""), ...)
    cat(paste(c("",
                "Results per cluster", 
                "-------------------", ""), 
              collapse = '\n'))
    print(as.data.frame(o$cluster, row.names = ""), ...)
    cat(paste(c("",
                "Non-redundancy per pair of clusters",
                "-----------------------------------", ""), 
              collapse = '\n'))
    .unquote(o$nonredundancy.clusterpairs, ...)
    cat(paste(c("",
                "RV per pair of clusters",
                "-----------------------", ""), 
              collapse = '\n'))
    .unquote(o$rv.clusterpairs, ...)
    # cat(paste(c("",
    #             "Salton's cosine per pair of clusters",
    #             "------------------------------------", ""), 
    #           collapse = '\n'))
    # .unquote(o$salton.clusterpairs, ...)
    o <- NULL
  }
  invisible(res)
}


