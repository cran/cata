#' Calculate within-cluster homogeneity
#'
#' Within a group of \code{N} consumers, the Homogeneity index lies between
#' \code{1/N} (no homogeneity) to \code{1} (perfect homogeneity). 
#' @name homogeneity
#' @aliases homogeneity
#' @usage homogeneity(X)
#' @param X three-way array; the \code{I, J, M} array has \code{I}
#' assessors, \code{J} products, \code{M} attributes where CATA data have values 
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
#' assessors, \code{J} products, \code{M} attributes where CATA data have values 
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
#' evaluateClusterQuality(bread$cata[1:8,,1:5], M = rep(1:2, each = 4))
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

#' Adjusted Rand index
#'
#' Calculate the adjusted Rand index between two sets of cluster memberships. 
#' @name ARI
#' @aliases ARI
#' @usage ARI(x, y, signif = FALSE, n = 1000)
#' @param x vector of cluster memberships (integers)
#' @param y vector of cluster memberships (integers)
#' @param signif conduct significance test; default is \code{FALSE}
#' @param n number of replicates in Monte Carlo significance test
#' @return ari adjusted Rand index
#' @return nari normalized adjusted Rand index 
#' @return sim.mean average value of null distribution (should be closed to zero)
#' @return sim.var variance of null distribution
#' @return pvalue P value of observed ARI (or NARI) value
#' @export
#' @encoding UTF-8
#' @references Hubert, L., & Arabie, P. (1985). Comparing partitions. 
#' \emph{Journal of Classification}, 2, 193–218.
#' \doi{10.1007/BF01908075}.
#' @references Qannari, E.M., Courcoux, P., & Faye, P. (2014). 
#' Significance test of the adjusted Rand index. Application to the free sorting 
#' task. \emph{Food Quality and Preference}, 32, 93-97. 
#' \doi{10.1016/j.foodqual.2013.05.005}.
#' 
#' @examples
#' x <- sample(1:3, 20, replace = TRUE)
#' y <- sample(1:3, 20, replace = TRUE)
#' 
#' ARI(x, y, signif = FALSE)
ARI <- function(x, y, signif = FALSE, n = 1000){
  getARI <- function(x, y){
    ux <- unique(x)
    uy <- unique(y)
    ub <- unique(ux, uy)
    M <- matrix(0, nrow=length(ub), ncol=length(ub))
    for(ix in 1:nrow(M)){
      for(iy in 1:ncol(M)){
        M[ix, iy] <- sum(x == ub[ix] & y == ub[iy])
      }
    }
    MM <- list(cellsums = c(M),
               rs = rowSums(M),
               cs = colSums(M),
               sum.ij = 0,
               sum.i = 0,
               sum.j = 0,
               cn = choose(length(x), 2))
    for(i in 1:length(MM$cellsums)){
      MM$sum.ij <- MM$sum.ij + choose(MM$cellsums[i], 2)
    }
    for(i in 1:length(MM$rs)){
      MM$sum.i <- MM$sum.i + choose(MM$rs[i], 2)
      MM$sum.j <- MM$sum.j + choose(MM$cs[i], 2)
    }
    
    ari <- (MM$sum.ij - (MM$sum.i * MM$sum.j / MM$cn)) /
      ((.5 * (MM$sum.i + MM$sum.j)) - (MM$sum.i * MM$sum.j / MM$cn))
    return(ari)
  }
  if(length(y) != length(x)){
    return("x and y must be equal lengths")
  }
  
  ari <- getARI(x, y)
  
  if(signif){
    xshuff <- yshuff <- matrix(0, nrow=n, ncol=length(x))
    ari.sim <- rep(0, n)
    for(i in 1:n){
      xshuff[i,] <- sample(x, length(x), replace=FALSE)
      yshuff[i,] <- sample(y, length(x), replace=FALSE)
      ari.sim[i] <- getARI(xshuff[i,], yshuff[i,])
    }
    nari.sim.mean <- mean(ari.sim)
    nari.sim.var <- stats::var(ari.sim)
    
    nari.sim <- (ari.sim - nari.sim.mean) / sqrt(nari.sim.var)
    nari <- (ari - nari.sim.mean) / sqrt(nari.sim.var)
    pval <- (sum(nari.sim > nari) + 1)/(n + 1)
    #pval.one <- (sum(ifelse(1-abc$ari.sim > 1, 1, abc$ari.sim) < ari))/(n + 1)
    return(list(ari = ari, 
                nari = nari, 
                sim.mean = nari.sim.mean,
                sim.var = nari.sim.var,
                pvalue = pval))
  } else {
    return(list(ari = ari))
  }
}




#' Plot variation in retained sensory differentiation
#'
#' Plot variation in retained sensory differentiation of cluster memberships
#' obtained from b-cluster analysis. This plot can be used to help the decision
#' of how many clusters to retain. 
#' @name selectionPlot
#' @aliases selectionPlot
#' @usage selectionPlot(x, pctB = NULL, x.input = "deltaB", indx = NULL, 
#' ylab = "change in B (K to G)", xlab = NULL)
#' @param x input vector which is either deltaB (default; change 
#' in sensory differentiation retained) or B (sensory differentiation 
#' retained) if \code{x.input} is \code{"B"}
#' @param pctB vector of percentage of the total sensory differentiation retained
#' @param x.input indicates what \code{x} is; either \code{"deltaB"} (default) 
#' or \code{B}.
#' @param indx numeric value indicating which point(s) to emphasize
#' @param ylab label shown on y axis and at selection point
#' @param xlab label for points along x axis
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Meyners, M., Varela, P., & Næs, T. (2022). 
#' Clustering consumers based on product discrimination in check-all-that-apply 
#' (CATA) data. \emph{Food Quality and Preference}, 104564. 
#' \doi{10.1016/j.foodqual.2022.104564}.
#'  
#' @examples
#' set.seed(123)
#' G2 <- bcluster.n(bread$cata[1:8, , 1:5], G = 2, runs = 3)
#' G3 <- bcluster.n(bread$cata[1:8, , 1:5], G = 3, runs = 3)
#' G4 <- bcluster.n(bread$cata[1:8, , 1:5], G = 4, runs = 3)
#' 
#' best.indx <- c(which.max(unlist(lapply(G2, function(x) x$retainedB))),
#'                which.max(unlist(lapply(G3, function(x) x$retainedB))),
#'                which.max(unlist(lapply(G4, function(x) x$retainedB))))
#'                
#' G1.bc <- barray(bread$cata[1:8, , 1:5])
#' G1.B <- getb(G1.bc[,,1,], G1.bc[,,2,])
#' BpctB <- data.frame(retainedB = c(G1.B, 
#'                                   G2[[best.indx[1]]]$retainedB, 
#'                                   G3[[best.indx[2]]]$retainedB,
#'                                   G4[[best.indx[3]]]$retainedB))
#' BpctB$pctB <- 100*BpctB$retainedB / G2[[1]]$totalB
#' BpctB$deltaB <- 
#'            c(100*(1-BpctB$retainedB[-nrow(BpctB)] / BpctB$retainedB[-1]), NA)
#' BpctB <- BpctB[-nrow(BpctB),]
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mar = rep(5,4))
#' selectionPlot(BpctB$deltaB, BpctB$pctB, indx = 2)
#' par(opar)
selectionPlot <- function(x, pctB = NULL, x.input = "deltaB", indx = NULL,
                          ylab = "change in B (K to G)", xlab = NULL){
  #requireNamespace("latex2exp", quietly = TRUE)
  x <- as.numeric(x)
  deltaB <- x # by default, x is deltaB
  if(x.input %in% "B"){
    deltaB <- 100*(1-x[-length(x)]/x[-1])
  }
  pctB <- as.numeric(pctB)
  nG <- length(deltaB)
  ylim.deltaB <- c(0, round(max(deltaB)+5,-1))
  ylim.pctB <- c(0, round(max(pctB)+5,-1))
  xlim <- c(.5, length(deltaB)+.5)
  plot(1:nG, deltaB, ylim = ylim.deltaB,
       xlab = "Change in number of clusters", 
       ylab = ylab, #TeX("$\\Delta B_{K \\rightarrow G}$%")
       type = "n", axes = FALSE) #
  if(is.null(xlab)){
    xlab <- paste0(2:(nG+1), "->", 1:nG)
  }
  graphics::axis(1, labels = xlab, # TeX(paste0(2:(nG+1), "$$\\rightarrow$$", 1:nG))
       at = 1:nG, las = 1, lwd = 4)
  graphics::axis(2, at = seq(0, ylim.deltaB[2], by = 10), 
       labels = seq(0, ylim.deltaB[2], by = 10), lwd = 4, las = 1)
  graphics::lines(x = 1:nG, y = deltaB, type = "b", pch = 16, lwd = 4)
  graphics::axis(4, labels = seq(0, ylim.pctB[2], by = 10), 
       at= seq(0, ylim.pctB[2], by = 10)*ylim.deltaB[2]/ylim.pctB[2], las=1)
  graphics::axis(3, labels = 1:nG, at = 1:nG, las = 1)
  graphics::mtext("Retained sensory differentiation (%B)", side = 4, line = 3)
  graphics::mtext("Number of clusters", side = 3, line = 3)
  graphics::lines(x = 1:nG, y = pctB*ylim.deltaB[2]/ylim.pctB[2], type = "b", pch = 17) # lty = "dashed", 
  if(!is.null(indx[1])){
    indx <- as.numeric(indx)
    for (i in seq_along(indx)){
      graphics::points(x=indx[i], y=deltaB[indx[i]], cex = 2.5, pch = 1)
      graphics::points(x=indx[i], y=pctB[indx[i]]*ylim.deltaB[2]/ylim.pctB[2], cex = 2.5, pch = 1)
      graphics::text(x = indx[i] + .1, y = mean(deltaB[indx[i]]) + .5, 
           labels = ylab, #TeX("$\\Delta B_{K \\rightarrow G}$%"), 
           pos = 4)
      graphics::text(x = indx[i] + .1, y = (mean(pctB[indx[i]]) - 1)*ylim.deltaB[2]/ylim.pctB[2], 
           labels = "%B", pos = 4)
    }
  }
  invisible()
}

