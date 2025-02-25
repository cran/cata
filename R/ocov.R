#' Permutation tests for CATA data
#'
#' Permutation tests for check-all-that-apply (CATA) data following the 
#' 'one citation, one vote' principle. Returns CATA frequency and percentage tables per 
#' condition and permutation test results specified.
#'
#' @name madperm
#' @usage madperm(X, B = 99, seed = .Random.seed, tests = 1:5, alpha = 0.05, control.fdr = FALSE, 
#' verbose = FALSE)
#' @param X a three-way (or four-way) array with \eqn{I} assessors, 
#' (\eqn{R} conditions/timepoints,) \eqn{P} products, \eqn{T} terms with data valued 
#' \code{0} (not checked) and \code{1} (checked)
#' @param B permutations in null distribution; ensure \code{B} is
#' larger than the default value (\code{99}) when conducting real analyses
#' @param seed specify a numeric seed for reproducibility; if not provided, a random 
#' seed is generated 
#' @param tests numeric vector specifying which tests to conduct; default 
#' \code{1:5} (all tests): (\code{1}) multivariate (global) test; (\code{2}) univariate tests; 
#' (\code{3}) elementwise tests; (\code{4}) pairwise multivariate tests; (\code{5}) 
#' pairwise univariate tests
#' @param alpha Type I error rate (default: \eqn{\alpha} = \code{0.05}) 
#' @param control.fdr control False Discovery Rate (using Benjamini-Hochberg (BH) step-up 
#' procedure)? (default: \code{FALSE})
#' @param verbose return null distribution(s) and function call? (default: \code{FALSE})
#' @return list, one per condition:
#' \itemize{
#' \item{\code{CATA.table} : table of CATA citation percentages (\eqn{P \times T})}
#' \item{\code{CATA.freq} : CATA frequency table (\eqn{P \times T})}
#' \item{Permutation test results specified by the \code{tests} parameter
#' \enumerate{
#' \item{\code{Global.Results} : list of multivariate (global) results}
#' \item{\code{Univariate.Results} : list of \eqn{T} univariate results}
#' \item{\code{Elementwise.Results} : list of \eqn{PT} elementwise results}
#' \item{\code{Multivariate.Paired.Results} : list of \eqn{P(P-1)/2} multivariate paired results}
#' \item{\code{Univariate.Paired.Results} : list of \eqn{P(P-1)T/2} univariate paired results} }}}
#' also, if \code{verbose} is \code{TRUE}:
#' \itemize{
#' \item{\code{Null.Dist} list of null distributions for tests specified}
#' \item{\code{Call} : \code{madperm} function call}}
#' @export
#' @encoding UTF-8
#' @author J.C. Castura
#' @references
#' Chaya, C., Castura, J.C., & Greenacre, M.J. (2025). One citation, one vote! 
#' A new approach for analyzing check-all-that-apply (CATA) data in sensometrics, 
#' using L1 norm methods. \doi{doi:10.48550/arXiv.2502.15945}
#' @examples
#' data(bread)
#' # add product names
#' X <- bread$cata[1:100,,1:5]
#' dimnames(X)[[2]] <- paste0("P", dimnames(X)[[2]])
#' # permutation tests for the first 100 consumers and 5 attributes
#' # will be run with default parameter values for illustrative purposes only
#' res <- madperm(X, B = 99, seed = 123)
#' print(res) # inspect results
madperm <- function(X, B = 99, seed = .Random.seed, tests = 1:5, 
                    alpha = 0.05, control.fdr = FALSE, 
                    verbose = FALSE){
  if((length(tests)>5) || !all(tests %in% 1:5)){
    return(print("tests must all be in 1 to 5"))
  }
  if(verbose){
    # report times per segment
    systimes <- rep(NA, 8)
    systimes.indx <- 1
    systimes[systimes.indx] <- Sys.time()
    systimes.indx <- systimes.indx + 1
  }
  .madslice <- function(Xm, addNames = FALSE){
    madfast <- apply(Xm, 2, stats::mad, constant = 1)
    if(isTRUE(addNames)){
      o <- c(stats::median(madfast), madfast)
      if(is.null(names(Xm[1]))){
        names(o) <- c("Global", seq_along(Xm)) 
      } else { 
        names(o) <- c("Global", names(Xm)) 
      }
      return(o)
    } else {
      return(c(stats::median(madfast), madfast))
    }
  }
  .nameDims <- function(A){
    for(i in 1:length(dim(A))){
      if(all(is.null(dimnames(A)[[i]]))){
        dimnames(A)[[i]] <- 1:(dim(A)[i])
      }
    }
    return(A)
  }
  .p.ge <- function(y){ 
    mean(abs(y) >= abs(y[1]))
  }
  .prop2pct <- function(x, digits = 1){
    return(round(x*100, digits = digits))
  }
  .sweepMedians <- function(tbl, Margins=NULL){
    sweep(tbl, Margins, .marginMedians(tbl, Margins))
  }
  .marginMedians <- function(tbl, Margins = NULL){
    apply(tbl, Margins, stats::median, drop = FALSE)
  }
  .fdrResults <- function(p.value, alpha = 0.05, digits = 8){
    BH.value = alpha*(1:length(p.value))/length(p.value)
    Signif = p.value < BH.value
    if(sum(Signif)>0){
      # match p-values in both directions to account for ties
      last.indx <- max(which.max(Signif), # forward direction
                       length(Signif) - which.max(rev(Signif)) + 1) # backward direction
      Signif[1:last.indx] <- TRUE
    }
    return(list(p.value = round(p.value, digits = digits),
                BH.value = round(BH.value, digits = digits),
                Signif = Signif))
  }
  .array_num2int <- function(x){
    array(as.integer(x), dim=dim(x), dimnames=dimnames(x)) 
  }
  .vec_num2int <- function(x){
    stats::setNames(as.integer(x), names(x))
  }
  In <- Rn <- Jn <- Mn <- 0
  if(length(dim(X)) == 3){
    .array_3to4 <- function(X){
      newdim <- c(1,dim(X))[c(2,1,3,4)]
      newdimnames <- append(dimnames(X), list(NULL), after = 1)
      return(array(as.integer(X), dim = newdim, dimnames = newdimnames) )
    }
    X <- .array_3to4(X)
  }
  if(length(dim(X)) == 4){
    X <- .array_num2int(X)
    In <- dim(X)[1]
    Rn <- dim(X)[2]
    Jn <- dim(X)[3]
    Mn <- dim(X)[4]
    X <- .nameDims(X)
    Rnames <- dimnames(X)[[2]]
    Jnames <- dimnames(X)[[3]]
    JJgrid <- utils::combn(Jn,2)
    JJn <- ncol(JJgrid)
    JJnames.grid <- matrix(Jnames[JJgrid], ncol=ncol(JJgrid))
    JJnames <- paste0(JJnames.grid[1,], "-", JJnames.grid[2,])
    JJcoded <- paste0(JJgrid[1,], "-", JJgrid[2,])
    Mnames <- dimnames(X)[[4]]
  }
  if(Jn == 0){
    return("Input error: cannot calculate results without a CATA data array")
  }
  if(verbose){
    systimes[systimes.indx] <- Sys.time()
    systimes.indx <- systimes.indx + 1
  }
  PermRes <- apply(X, 2, function(Y){ list(
    CATA.table = apply(Y, 2:3, mean, drop = FALSE),
    CATA.freq = apply(Y, 2:3, sum)
  ) }, 
                   simplify=FALSE)
  names(PermRes) <- Rnames
  CATA.freq <- apply(X, 2:4, sum)[,,, drop=FALSE]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Permutation tests
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Part 1.Observed results
  if(verbose){
    systimes[systimes.indx] <- Sys.time()
    systimes.indx <- systimes.indx + 1
  }
  nulldist <- list()
  if(any(1:2 %in% tests)){
    nulldist$variatewise  <- array(NA, c(Rn, 1+B, 1+Mn), dimnames = list(
      paste0("rep_", 1:Rn),
      c("Observed", paste0("b", 1:B)),
      c("Global", Mnames)))
    nulldist$variatewise[,1,] <- t(apply(CATA.freq, #CATA.table, 
                                         1, .madslice)) # Rn x Mn
  }
  if(any(3 %in% tests)){
    nulldist$elementwise <- array(NA, c(Rn, 1+B, Jn, Mn), dimnames = list(
      paste0("rep_", 1:Rn),
      c("Observed", paste0("b", 1:B)),
      1:Jn,
      Mnames))
    nulldist$elementwise[,1,,] <- 
      .sweepMedians(CATA.freq, #CATA.table, 
                    Margins = c(1,3))
  }
  if(any(4:5 %in% tests)){
    nulldist$pairwise <- 
      array(NA, c(Rn, 1+B, 1+Mn, Jn*(Jn-1)/2), dimnames = list(
        paste0("rep_", 1:Rn),
        c("Observed", paste0("b", 1:B)), 
        c("Global", Mnames), 
        JJnames))
    for(jj in 1:JJn){
      nulldist$pairwise[,1,-1,jj] <- 
        t(apply(CATA.freq[ ,c(JJgrid[2,jj], JJgrid[1,jj]),, drop=FALSE], 1, diff))
    }
    nulldist$pairwise[,1,1,] <- 
      .marginMedians(abs(nulldist$pairwise[,1,-1,, drop=FALSE]), Margins = c(1,4))
  }
  # Part 2. Permutation results
  if(verbose){
    systimes[systimes.indx] <- Sys.time()
    systimes.indx <- systimes.indx + 1
  }
  Xb <- X # setup Xb once, here before starting the loop
  set.seed(seed)
  for(b in 1:B){
    # permute products within assessors
    Jorder <- replicate(In, sample(Jn, Jn, replace=FALSE))
    for(i in 1:In){
      Xb[i,,,] <- X[i,,Jorder[,i],, drop=FALSE] # product order is the same for all reps
    }
    # if(DEBUG && b==1){
    #   Xb <- X
    # }
    Xb.freq <- apply(Xb, 2:4, sum) 
    if(any(1:2 %in% tests)){
      nulldist$variatewise[,1+b,] <- 
        t(apply(Xb.freq, 1, .madslice)) # results: rep x attribute
      # check: 
      # t(apply(abs(sweep(drop(Xb.freq),2,apply(drop(Xb.freq),2,median))),2,median)) - 
      #   nulldist$variatewise[,1+b,][-1]
    }
    if(any(3 %in% tests)){
      nulldist$elementwise[,1+b,,] <- 
        .sweepMedians(Xb.freq, Margins = c(1,3))
    }
    if(any(4:5 %in% tests)){
      for(jj in 1:JJn){
        nulldist$pairwise[,1+b,-1,jj] <- 
          apply(Xb.freq[,c(JJgrid[2,jj], JJgrid[1,jj]),, drop=FALSE], c(1,3), 
                diff) # 1st - 2nd
      }
      nulldist$pairwise[,1+b,1,] <- 
        .marginMedians(abs(nulldist$pairwise[,1+b,-1,, drop=FALSE]), Margins = c(1,4))
    }
  }
  # Part 3. Calculate P values
  if(verbose){
    systimes[systimes.indx] <- Sys.time()
    systimes.indx <- systimes.indx + 1
  }
  pval <- list()
  if(any(1:2 %in% tests)){
    # P values for variatewise tests
    pval$global <- apply(nulldist$variatewise, c(1,3), .p.ge)
  }
  if(any(3 %in% tests)){
    pval$elementwise <- 
      apply(nulldist$elementwise, c(1,3,4), .p.ge)
  }
  if(any(4:5 %in% tests)){
    # P values for pairwise tests
    # Simplify to return matrix/matrix, but maybe will it also drop in the case of 1-d?
    pval$pairwise <- apply(nulldist$pairwise, c(1,3,4), .p.ge, simplify = TRUE)
  }
  # Summarize results
  if(verbose){
    systimes[systimes.indx] <- Sys.time()
    systimes.indx <- systimes.indx + 1
  }
  p.dig <- ceiling(max(1, log10(B+1))) # digits for p, BH values
  if(1 %in% tests){
    for(rindx in 1:Rn) {
      PermRes[[rindx]]$Global.Results <- 
        data.frame(
        MAD = nulldist$variatewise[rindx,1,1, drop=TRUE],
        MAD.pct = .prop2pct(nulldist$variatewise[rindx,1,1, drop=TRUE]/In),
        p.value = round(pval$global[rindx, 1], digits = p.dig),
        Signif = pval$global[rindx, 1] <= alpha )
    }
  }
  if(2 %in% tests){
    if(!control.fdr){
      for (rindx in 1:Rn){
        PermRes[[rindx]]$Univariate.Results <- 
          data.frame(
            Term.Indx = 1:Mn,
            Term = Mnames,
            MAD = nulldist$variatewise[rindx,1,-1, drop=TRUE],
            MAD.pct = .prop2pct(nulldist$variatewise[rindx,1,-1, drop=TRUE]/In),
            p.value = round(pval$global[rindx, -1], digits = p.dig), 
            Signif = pval$global[rindx, -1] <= alpha)
        PermRes[[rindx]]$Univariate.Results <- 
          PermRes[[rindx]]$Univariate.Results[
            order(pval$global[rindx, -1], decreasing = FALSE),]
      }
    } else {
      for (rindx in 1:Rn){
        PermRes[[rindx]]$Univariate.Results <- 
          data.frame(
            Term.Indx = 1:Mn,
            Term = Mnames,
            MAD = nulldist$variatewise[rindx,1,-1, drop=TRUE],
            MAD.pct = .prop2pct(nulldist$variatewise[rindx,1,-1, drop=TRUE]/In),
            p.value = pval$global[rindx, -1])
        PermRes[[rindx]]$Univariate.Results <- 
          PermRes[[rindx]]$Univariate.Results[
            order(PermRes[[rindx]]$Univariate.Results$p.value, decreasing = FALSE),]
        PermRes[[rindx]]$Univariate.Results <- 
            cbind(PermRes[[rindx]]$Univariate.Results[, -5], # replace p values in column 5
                  .fdrResults(PermRes[[rindx]]$Univariate.Results$p.value, digits = p.dig))
      }
    }
  }
  if(3 %in% tests){
    elements <- expand.grid(Product.Indx = 1:Jn, Term.Indx = 1:Mn)
    elements$Product <- Jnames[elements$Product.Indx]
    elements$Term <- Mnames[elements$Term.Indx]
    if(!control.fdr){
      for (rindx in 1:Rn){
        PermRes[[rindx]]$Elementwise.Results <- 
          data.frame(
            elements,
            MAD = c(nulldist$elementwise[rindx,1,,, drop=FALSE]),
            MAD.pct = .prop2pct(c(nulldist$elementwise[rindx,1,,, drop=FALSE])/In),
            p.value = round(c(pval$elementwise[rindx,,]), digits = p.dig),
            Signif = c(pval$elementwise[rindx,,]) <= alpha)
        PermRes[[rindx]]$Elementwise.Results <- 
          PermRes[[rindx]]$Elementwise.Results[
            order(c(pval$elementwise[rindx,,]), decreasing = FALSE),]
      }
    } else {
      for (rindx in 1:Rn){
        PermRes[[rindx]]$Elementwise.Results <- 
          data.frame(
            elements,
            MAD = c(nulldist$elementwise[rindx,1,,, drop=FALSE]),
            MAD.pct = .prop2pct(c(nulldist$elementwise[rindx,1,,, drop=FALSE])/In),
            p.value = c(pval$elementwise[rindx,,])) # p values are not rounded yet
        PermRes[[rindx]]$Elementwise.Results <- 
          PermRes[[rindx]]$Elementwise.Results[
            order(c(pval$elementwise[rindx,,]), decreasing = FALSE),]
        PermRes[[rindx]]$Elementwise.Results <- 
          cbind(PermRes[[rindx]]$Elementwise.Results[, -7],  # replace p values in column 7
                .fdrResults(PermRes[[rindx]]$Elementwise.Results$p.value, digits = p.dig))
      }
    }
  }
  if(4 %in% tests){
    if(!control.fdr){
    for (rindx in 1:Rn){
      PermRes[[rindx]]$Multivariate.Paired.Results <- 
        data.frame(
          Pair.Indx = 1:JJn,
          Pair = JJcoded,
          Pair.name = JJnames,
          MAD = as.vector(t(abs(nulldist$pairwise[rindx,1,1,]))),
          MAD.pct = .prop2pct(as.vector(t(abs(nulldist$pairwise[rindx,1,1,])))/In),
          p.value = round(unlist(pval$pairwise[rindx,1,]), digits = p.dig),
          Signif = unlist(pval$pairwise[rindx,1,]) <= alpha)
      PermRes[[rindx]]$Multivariate.Paired.Results <- 
        PermRes[[rindx]]$Multivariate.Paired.Results[
          order(PermRes[[rindx]]$Multivariate.Paired.Results$p.value, 
                decreasing = FALSE),]
    }
      } else {
        PermRes[[rindx]]$Multivariate.Paired.Results <- 
          data.frame(
            Pair.Indx = 1:JJn,
            Pair = JJcoded,
            Pair.name = JJnames,
            MAD = as.vector(t(abs(nulldist$pairwise[rindx,1,1,]))),
            MAD.pct = .prop2pct(as.vector(t(abs(nulldist$pairwise[rindx,1,1,])))/In),
            p.value = unlist(pval$pairwise[rindx,1,]))
        PermRes[[rindx]]$Multivariate.Paired.Results <- 
          PermRes[[rindx]]$Multivariate.Paired.Results[
            order(PermRes[[rindx]]$Multivariate.Paired.Results$p.value, 
                  decreasing = FALSE),]
        PermRes[[rindx]]$Multivariate.Paired.Results <- 
          cbind(PermRes[[rindx]]$Multivariate.Paired.Results[, -6],  # replace p values in column 6
                .fdrResults(PermRes[[rindx]]$Multivariate.Paired.Results$p.value, digits = p.dig))
    }
  }
  if(5 %in% tests){
    pair.elements <- expand.grid(Pair.Indx = 1:JJn, Term.Indx = 1:Mn)
    pair.elements$Pair <- JJcoded[pair.elements$Pair.Indx]
    pair.elements$Pair.names <- JJnames[pair.elements$Pair.Indx]
    pair.elements$Term <- Mnames[pair.elements$Term.Indx]
    if(!control.fdr){
      for (rindx in 1:Rn){
          PermRes[[rindx]]$Univariate.Paired.Results <- 
        data.frame(
          pair.elements,
          MAD = c(t(nulldist$pairwise[rindx,1,-1,])),
          MAD.pct = .prop2pct(c(t(nulldist$pairwise[rindx,1,-1,]))/In),
          p.value = round(c(t(pval$pairwise[rindx,-1,])), digits = p.dig), 
          Signif = c(t(pval$pairwise[rindx,-1,])) <= alpha)
      PermRes[[rindx]]$Univariate.Paired.Results <- 
        PermRes[[rindx]]$Univariate.Paired.Results[
          order(c(t(pval$pairwise[rindx,-1,])), decreasing = FALSE),]
      }
    } else {
      for (rindx in 1:Rn){
        PermRes[[rindx]]$Univariate.Paired.Results <- 
          data.frame(
            pair.elements,
            MAD = c(t(nulldist$pairwise[rindx,1,-1,])),
            MAD.pct = .prop2pct(c(t(nulldist$pairwise[rindx,1,-1,]))/In),
            p.value = c(t(pval$pairwise[rindx,-1,])))
        PermRes[[rindx]]$Univariate.Paired.Results <- 
          PermRes[[rindx]]$Univariate.Paired.Results[
            order(c(t(pval$pairwise[rindx,-1,])), decreasing = FALSE),]
        PermRes[[rindx]]$Univariate.Paired.Results <- 
          PermRes[[rindx]]$Univariate.Paired.Results[
            order(PermRes[[rindx]]$Univariate.Paired.Results$p.value, decreasing = FALSE),]
        PermRes[[rindx]]$Univariate.Paired.Results <- 
          cbind(PermRes[[rindx]]$Univariate.Paired.Results[, -8], # replace p values in column 8
                .fdrResults(PermRes[[rindx]]$Univariate.Paired.Results$p.value, digits = p.dig))
      }
    }
  }
  if(verbose){
    systimes[systimes.indx] <- Sys.time()
    systimes.indx <- systimes.indx + 1
  }
  if(Rn == 1){
    o <- PermRes[[1]]
    o$CATA.table <- .prop2pct(o$CATA.table)
  } else {
    o <- list(Permutation.Tests = PermRes)
    for(rindx in 1:Rn){
      o$Permutation.Tests[[rindx]]$CATA.table <- 
        .prop2pct(o$Permutation.Tests[[rindx]]$CATA.table)
    }
  }
  if(isTRUE(verbose)){
    o$Call <- match.call()
    o$Null.Dist <- list()
    if(1 %in% tests) o$Null.Dist$Global <- nulldist$variatewise[,,1, drop=TRUE]
    if(2 %in% tests) o$Null.Dist$Univariate <- nulldist$variatewise[,,-1, drop=TRUE]
    if(3 %in% tests) o$Null.Dist$Elementwise <- drop(nulldist$elementwise)
    if(4 %in% tests) o$Null.Dist$Multivariate.Paired <- nulldist$pairwise[,,1,, drop=TRUE] 
    if(5 %in% tests) o$Null.Dist$Univariate.Paired <- nulldist$pairwise[,,-1,, drop=TRUE]
  }
  if(verbose){
    systimes[systimes.indx] <- Sys.time()
    o$Runtimes <- c(systimes[-1] - systimes[-length(systimes)])
    names(o$Runtimes) <- c("Init", "CATA.tables", "Observed", "Permutations",
                           "Pvalues", "Summarize", "Return")
  }
  return(o)
}

#' MAD distances between objects
#'
#' Computes and returns inter-object median of absolute deviations (MADs) based 
#' differences. 
#'
#' @name mad.dist
#' @usage mad.dist(X)
#' @param X objects-by-terms matrix
#' @return object of class \code{dist} giving inter-object MAD distances 
#' @export
#' @encoding UTF-8
#' @author J.C. Castura
#' @references 
#' Chaya, C., Castura, J.C., & Greenacre, M.J. (2025). One citation, one vote! 
#' A new approach for analyzing check-all-that-apply (CATA) data in sensometrics, 
#' using L1 norm methods. \doi{doi:10.48550/arXiv.2502.15945}
#' @examples
#' data(bread)
#' CATA.freq <- apply(bread$cata, 2:3, sum)
#' # median-center columns (attributes)
#' CATA.swept <- sweep(CATA.freq, 2, apply(CATA.freq, 2, median))
#' # cluster analysis of products using complete linkage
#' dist.Products <- mad.dist(CATA.swept)
#' plot(as.dendrogram(hclust(dist.Products, method = "complete")), main = "Product clusters")
#' # cluster analysis of attributes using complete linkage
#' dist.Att <- mad.dist(t(CATA.swept))
#' plot(as.dendrogram(hclust(dist.Att, method = "complete")), main = "Attribute clusters")
mad.dist <- function(X) {
  Xcomb <- utils::combn(nrow(X), 2)
  o <- matrix(0, nrow = nrow(X), ncol = nrow(X), dimnames = list(rownames(X), 
                                                                 rownames(X)))
  for (i in 1:ncol(Xcomb)) {
    o[Xcomb[1, i], Xcomb[2, i]] <- 
      o[Xcomb[2, i], Xcomb[1, i]] <- 
      stats::median(abs(apply(X[c(Xcomb[1, i], Xcomb[2, i]), ], 2, diff)))
  }
  return(stats::as.dist(o))
}

