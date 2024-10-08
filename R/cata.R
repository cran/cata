#' Wrapper function for b-cluster analysis
#'
#' By default, \code{bcluster} calls a function to perform b-cluster analysis 
#' by a non-hierarchical iterative ascent algorithm, then inspects results if 
#' there are multiple runs. 
#' 
#' @name bcluster
#' @aliases bcluster
#' @usage bcluster(X, inspect = TRUE, inspect.plot = TRUE, algorithm = "n", 
#' measure = "b", G = NULL, M = NULL, max.iter = 500, X.input = "data", 
#' tol = exp(-32), runs = 1, seed = 2021)
#' @param X three-way array with \code{I} assessors, \code{J} products, 
#' \code{M} attributes where CATA data have values \code{0} (not checked) and 
#' \code{1} (checked)
#' @param inspect default (\code{TRUE}) calls the \code{\link[cata]{inspect}}
#' function to evaluate all solutions (when \code{runs>1})
#' @param inspect.plot default (\code{TRUE}) plots results from the 
#' \code{\link[cata]{inspect}} function 
#' @param algorithm default is \code{n} for non-hierarchical; \code{h} for 
#' hierarchical 
#' @param measure default is \code{b} for the \code{b}-measure; \code{Q} for 
#' Cochran's Q test 
#' @param G number of clusters (required for non-hierarchical algorithm)
#' @param M initial cluster memberships
#' @param max.iter maximum number of iteration allowed (default \code{500})
#' @param X.input available only for non-hierarchical algorithm; its value is
#' either \code{"data"} (default) or \code{"bc"} if \code{X} is
#' obtained from the function \code{\link[cata]{barray}} 
#' @param tol non-hierarchical algorithm stops if variance over 5 iterations is
#' less than \code{tol} (default: \code{exp(-32)})
#' @param runs number of runs (defaults to \code{1})
#' @param seed for reproducibility (default is \code{2021})
#' @return list with elements:
#' \itemize{
#' \item{\code{runs} : b-cluster analysis results from 
#' \code{\link[cata]{bcluster.n}} or \code{\link[cata]{bcluster.h}} 
#' (in a list if \code{runs>1})}
#' \item{\code{inspect} : result from \code{\link[cata]{inspect}} (the plot from 
#' this function is rendered if \code{inspect.plot} is \code{TRUE})}}
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Meyners, M., Varela, P., & Næs, T. (2022). 
#' Clustering consumers based on product discrimination in check-all-that-apply 
#' (CATA) data. \emph{Food Quality and Preference}, 104564. 
#' \doi{10.1016/j.foodqual.2022.104564}.
#' @examples
#' data(bread)
#' 
#' # b-cluster analysis on the first 8 consumers and the first 5 attributes
#' (b1 <- bcluster(bread$cata[1:8,,1:5], G=2, seed = 123))
#' # Since the seed is the same, the result will be identical to
#' # (b2 <- bcluster.n(bread$cata[1:8,,1:5], G=2, seed = 123))
#' b3 <- bcluster(bread$cata[1:8,,1:5], G=2, runs = 5, seed = 123)
bcluster <- function(X, inspect = TRUE, inspect.plot = TRUE,
                     algorithm = "n", measure = "b",
                     G = NULL, M = NULL, max.iter = 500, X.input = "data", 
                     tol = exp(-32), runs = 1, seed = 2021){
  if(is.null(X)) return(print("X must be an array"))
  if(algorithm %in% c("h", "1")){
    if(X.input %in% "bc"){
      return(invisible("X.input 'bc' not available for hierarchical algorithm"))
    }
    out <- bcluster.h(X = X, measure = measure, runs = runs, 
                      seed = seed)
  }
  if(algorithm %in% c("n", "2")){
    if(is.null(G)){
      return(print("G must be specified"))
    } 
    if(length(G) > 1){
      return(print("G cannot be longer than one"))
    }
    out <- bcluster.n(X = X, G = G, M = M, measure = measure, 
                      max.iter = max.iter, X.input = X.input,
                      tol = tol, runs = runs,
                      seed = seed)
  }
  out <- list(runs = out)
  if((runs > 1) & inspect){
    if(algorithm %in% "h"){
      if(is.null(G)){
        print("Number of groups (G) not provided. Defaulting to G=2 groups.")
        G <- 2
      }
    }
    print(out$inspect <- inspect(out$runs, G = G, inspect.plot = inspect.plot))
  }
  invisible(out)
}

#' b-cluster analysis by hierarchical agglomerative strategy
#'
#' Perform b-clustering using the hierarchical agglomerative clustering 
#' strategy. 
#' @name bcluster.h
#' @aliases bcluster.h
#' @usage bcluster.h(X, measure = "b", runs = 1, seed = 2021)
#' @param X three-way array; the \code{I, J, M} array has \code{I}
#' assessors, \code{J} products, \code{M} attributes where CATA data have values 
#' \code{0} (not checked) and \code{1} (checked)
#' @param measure currently only \code{b} (the \code{b}-measure) is implemented 
#' @param runs number of runs (defaults to \code{1}; use a higher number of
#' runs for a real application)
#' @param seed for reproducibility (default is \code{2021})
#' @return An object of class \code{hclust} from hierarchical b-cluster 
#' analysis results (a list of such objects if \code{runs>1}), where each \code{hclust} 
#' object has the structure described in \code{\link[stats]{hclust}} as well as 
#' the item \code{retainedB} (a vector indicating the retained sensory 
#' differentiation at each iteration (merger)).
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Meyners, M., Varela, P., & Næs, T. (2022). 
#' Clustering consumers based on product discrimination in check-all-that-apply 
#' (CATA) data. \emph{Food Quality and Preference}, 104564. 
#' \doi{10.1016/j.foodqual.2022.104564}.
#' @examples
#' data(bread)
#' 
#' # hierarchical b-cluster analysis on first 8 consumers and first 5 attributes
#' b <- bcluster.h(bread$cata[1:8,,1:5])
#' 
#' plot(as.dendrogram(b), 
#'   main = "Hierarchical b-cluster analysis", 
#'   sub = "8 bread consumers on 5 attributes")
bcluster.h <- function(X, measure = "b", runs = 1, seed = 2021){
  ## DEBUG
  # X = bread$cata[1:14,,1]
  # measure = "b"
  # runs = 1
  # seed = 2021
  verbose = FALSE
  bcluster.call <- match.call()
  if(measure != "b"){
    return(print("Only the b-measure is implemented"))
  } 
  ### functions
  # get.gmem gets group members
  #  g.indx is the index for the group
  #  curr.df current allocations in the first 2 columns
  get.gmem <- function(g.indx, curr.df){
    return(curr.df$start.g[curr.df$curr.g == g.indx])
  }
  # init.candidates
  #   cnd1.indx : all members IDs of group 1
  #   cnd2.indx : all members IDs of group 2
  #   X.bc      : bc array
  #   b         : vector of b-measures
  #   sensDiff  : mN columns with 1 if differentiated and 0 if not
  init.candidates <- function(cnd1.indx, cnd2.indx, X.bc, b, sensDiff, 
                              X.bc.dim = NULL){
    # DEBUG
    # cnd1.indx = df.cnd$cnd1
    # cnd2.indx = df.cnd$cnd2
    # X.bc = X.bc
    # b = df.curr$b
    # sensDiff = df.curr[,-c(1:3)]
    # X.bc.dim = X.bc.dim
    if(!is.null(X.bc.dim[1])){
      if(!all.equal(dim(X.bc), X.bc.dim)){
        X.bc <- array(X.bc, X.bc.dim)
      }  
    }
    new.b <- getb(X.bc[c(cnd1.indx, cnd2.indx),,1,], 
                  X.bc[c(cnd1.indx, cnd2.indx),,2,],
                  oneI = ifelse(length(c(cnd1.indx, cnd2.indx))==1, TRUE, FALSE),
                  oneM = ifelse(X.bc.dim[length(X.bc.dim)]==1, TRUE, FALSE))
    if(X.bc.dim[length(X.bc.dim)]==1){
      sensDiff <- matrix(sensDiff, ncol=1)
      sensDiff.o <- sum(colSums(matrix(sensDiff[c(cnd1.indx, cnd2.indx), ], 
                                       ncol=1))>1)
    } else {
      sensDiff.o <- sum(colSums(sensDiff[c(cnd1.indx, cnd2.indx), ])>1)
    }
    return(c(new.b, new.b - b[cnd1.indx] - b[cnd2.indx],
             sensDiff.o))
  }
  # check if id being merged is a singleton (or not)
  isSingleton <- function(id, g.curr){
    return(sum(g.curr$curr.g %in% id)==1)
  }
  # update C and # differentiating attributes
  # update.candidates
  #   cnd1.indx : keys corresponding to all members of group 1
  #   cnd2.indx : keys corresponding to all members of group 2
  #   X.bc      : bc array
  #   b         : vector of b-measures
  #   sensDiff  : mN columns with 1 if differentiated and 0 if not
  update.candidates <- function(cnd1.indx, cnd2.indx, X.bc, 
                                group.keys, b, sensDiff, X.bc.dim = NULL){
    if(!is.null(X.bc.dim[1])){
      if(!all.equal(dim(X.bc), X.bc.dim)){
          X.bc <- array(X.bc, X.bc.dim)
      }  
    }
    g1.keys <- group.keys$start.g[group.keys$curr.g == cnd1.indx]
    g2.keys <- group.keys$start.g[group.keys$curr.g == cnd2.indx]
    new.b <- getb(X.bc[c(g1.keys, g2.keys),,1,], 
                  X.bc[c(g1.keys, g2.keys),,2,],
                  oneI = ifelse(length(c(g1.keys, g2.keys))==1, TRUE, FALSE),
                  oneM = ifelse(X.bc.dim[length(X.bc.dim)]==1, TRUE, FALSE))
    if(X.bc.dim[length(X.bc.dim)]==1){
      sensDiff <- matrix(sensDiff, ncol=1)
      sensDiff.o <- sum(colSums(matrix(sensDiff[c(group.keys$curr.g 
                                           %in% c(cnd1.indx, cnd2.indx)), ],
                                       ncol=1))>1)
    } else {
      sensDiff.o <- sum(colSums(sensDiff[c(group.keys$curr.g 
                                           %in% c(cnd1.indx, cnd2.indx)), ])>1)
    }
    return(c(new.b, 
             new.b - b[group.keys$start.g == cnd1.indx] - 
               b[group.keys$start.g == cnd2.indx],
             sensDiff.o))
  }
  # create hclust merge matrix from M, which is my progress.track object
  write.merge <- function(M){
    out <- M
    last.use <- rep(NA, nrow(M))
    for(i in 1:nrow(M)){
      if(!is.na(last.use[abs(M[i,1])])){
        out[i,1] <- last.use[M[i,1]]
      }
      if(!is.na(last.use[abs(M[i,2])])){
        out[i,2] <- last.use[M[i,2]]
      }
      last.use[abs(M[i,1])] <- i
      last.use[abs(M[i,2])] <- i
    }
    out <- as.matrix(out)
    colnames(out) <- NULL
    return(out)
  }
  ## end of functions
  
  nI <- dim(X)[1]
  nJ <- dim(X)[2]
  if(length(dim(X))==2){
    X <- array(X, c(dim(X)[1], dim(X)[2], 1),
               dimnames = list(dimnames(X)[[1]], dimnames(X)[[2]], "data"))
  } 
  nM <- dim(X)[3]
  X.bc <- barray(X)
  if(length(dim(X.bc))==3){
    X.bc <- array(X.bc, c(dim(X.bc)[1], dim(X.bc)[2], 2, 1),
               dimnames = list(dimnames(X.bc)[[1]], dimnames(X.bc)[[2]],
                              letters[2:3], "data"))
  } 
  nJJ <- dim(X.bc)[2]
  
  if(is.null(dimnames(X.bc)[[1]])) dimnames(X.bc)[[1]] <- 1:nI
  if(is.null(dimnames(X.bc)[[2]])) dimnames(X.bc)[[2]] <- 1:nJJ
  if(is.null(dimnames(X.bc)[[3]])) dimnames(X.bc)[[2]] <- c("n10", "n01")
  if(is.null(dimnames(X.bc)[[4]])) dimnames(X.bc)[[4]] <- 1:nM
  
  X.bc.dim <- dim(X.bc)
  
  progress.track <- data.frame(iter = 1:(nI-1), 
                               cons1 = NA, cons2 = NA, 
                               isSingleton1 = NA, isSingleton2 = NA,
                               C = NA,
                               n.tiebreak.sensDiff = NA, 
                               n.tiebreak.rand = NA)
  # df.curr
  #   $start.g   consumer start group indx
  #   $curr.g    current group indx
  #   $b         sensory differentiation of this group (NA if merged)
  #   $sensDif   an nM-column matrix for current group indicating that 1+ consumer 
  #              differentiates products on this attribute (NA if merged)
  df.curr <- data.frame(start.g = 1:nI, 
                        curr.g = 1:nI,
                        b = rowSums(abs(X.bc[,,1,]-X.bc[,,2,])),
                        sensDiff = ((apply(X.bc, c(1,4), sum) %% nJJ != 0)*1))
  totB <- sum(df.curr$b)
  
  # gMat
  # $cnd1    candidate group 1 for possible merger       
  # $cnd2    candidate group 2 for possible merger       
  # $C       C, sensory differentiation lost by merging candidate groups       
  # $sensDif number of attributes for which 1+ consumer differentiates products      
  df.cnd <- as.data.frame(cbind(t(utils::combn(nI, 2)), 
                                matrix(NA, ncol =3, nrow=nI*(nI-1)/2)))
  colnames(df.cnd) <- c("cnd1", "cnd2", "b", "C", "sensDif")
  
  df.cnd[, 3:5] <- t(mapply(init.candidates, 
                            cnd1.indx = df.cnd$cnd1, cnd2.indx = df.cnd$cnd2, 
                            MoreArgs = list(X.bc = X.bc, b = df.curr$b, 
                                            sensDiff = df.curr[,-c(1:3)],
                                            X.bc.dim = X.bc.dim)))
  if(runs>1){
    # cache calculated start values for future runs
    cache <- list(progress.track = progress.track,
                  df.cnd = df.cnd,
                  df.curr = df.curr)
  }
  hclust.list <- list()
  if(verbose){
    track.list <- list()
  }
  for(this.run in 1:runs){
    set.seed(seed + this.run)
    if(this.run > 1){
      # restore from cache
      progress.track <- cache$progress.track
      df.cnd <- cache$df.cnd
      df.curr <- cache$df.curr
    }
    # track output
    hclust.list[[this.run]] <- list()
    if(verbose){
      hclust.list[[this.run]]$track.list <-
        # track.list[[this.run]] <- list()
        # track.list[[this.run]][[1]] <- 
        list(df.cnd = df.cnd, df.curr = df.curr)
      hclust.list[[this.run]]$track.list.by.iter <- list()
    }
    for(iter in 1:(nI-1)){
      if(verbose){
        # Write state to track.list  
        hclust.list[[this.run]]$track.list.by.iter[[iter+1]] <- 
          # track.list[[this.run]][[iter+1]] <- 
          list(df.cnd = df.cnd, 
               df.curr = df.curr)
      }
      # Find groups that the maximize C
      best.indx <- which(df.cnd$C == max(df.cnd$C))
      if(length(best.indx)>1){
        # Tie-breaking required
        # 1. Maximize # attributes with sensory differentiation in both groups
        best.indx2 <- best.indx[which(df.cnd$sensDif[best.indx] == 
                                        max(df.cnd$sensDif[best.indx]))]
        # 2. Still ties to select at random
        g.indx <- as.numeric(sample(as.character(best.indx2), 1))
      } else {
        g.indx <- best.indx2 <- best.indx
      }
      # Update progress track
      progress.track[iter,] <- 
        c(iter, 
          df.cnd$cnd1[g.indx], df.cnd$cnd2[g.indx],
          isSingleton(df.cnd$cnd1[g.indx], df.curr[,1:2]),
          isSingleton(df.cnd$cnd2[g.indx], df.curr[,1:2]),
          df.cnd$C[g.indx], length(best.indx)-1, length(best.indx2)-1)
      # Post-merge updates
      # Update the current group of second group (progress.track$cons2[iter])
      df.curr$curr.g[which(df.curr$curr.g == progress.track$cons2[iter])] <-
        progress.track$cons1[iter]
      # Update the b-measure of new group
      df.curr$b[which(df.curr$curr.g == progress.track$cons1[iter])] <-
        df.cnd$b[g.indx]
      # Update number of attributes differentiated by new group
      df.curr[which(df.curr$start.g == 
                      progress.track$cons1[iter]), -c(1:3)] <-
        ifelse(nM > 1, (colSums(df.curr[
          which(df.curr$curr.g %in% 
                  c(progress.track$cons1[iter],
                    progress.track$cons2[iter])), -c(1:3)])>0)*1,
          (colSums(matrix(df.curr[
            which(df.curr$curr.g %in% 
                    c(progress.track$cons1[iter],
                      progress.track$cons2[iter])), -c(1:3)], ncol=nM))>0)*1)
      # Remove rows for second group from df.cnd
      df.cnd <- df.cnd[-c(which(df.cnd$cnd1 == progress.track$cons2[iter]),
                          which(df.cnd$cnd2 == progress.track$cons2[iter])),]
      # Update comparisons involving the new merged group
      update.indx <- sort(c(which(df.cnd$cnd1 == progress.track$cons1[iter]),
                            which(df.cnd$cnd2 == progress.track$cons1[iter])))
      df.cnd[update.indx, -c(1:2)] <-
        t(mapply(update.candidates, 
                 cnd1.indx = df.cnd$cnd1[update.indx],  
                 cnd2.indx = df.cnd$cnd2[update.indx],
                 MoreArgs = 
                   list(X.bc = X.bc, 
                        group.keys = df.curr[,1:2],
                        b = df.curr$b,
                        sensDiff = df.curr[, -c(1:3)],
                        X.bc.dim = X.bc.dim)))  
    }
    # All iterations complete so create the hclust object
    hc.obj <- structure(list(
      merge = write.merge(
        cbind(progress.track$cons1*c(1,-1)[progress.track$isSingleton1+1],
              progress.track$cons2*c(1,-1)[progress.track$isSingleton2+1])),
      height = -1 * progress.track$C,
      order = progress.track$iter,
      labels = as.character(dimnames(X.bc)[[1]]), 
      retainedB = totB + cumsum(progress.track$C),
      method = paste0("b-cluster analysis (", measure, "-measure)"), 
      call = bcluster.call, 
      seed = seed + this.run,
      dist.method = measure), 
      class = "hclust")
    hc.obj$order <- unname(stats::order.dendrogram(stats::as.dendrogram(hc.obj)))
    # if(!(stats::.validity.hclust(hc.obj))){
    #   print("Dendrogram failed")
    # }
    class(hc.obj) <- "hclust"
    hclust.list[[this.run]] <- hc.obj
    if(verbose){
      # Write state to track.list  
      hclust.list[[this.run]]$track.list.by.iter[[
        length(hclust.list[[this.run]]$track.list.by.iter)]]$progress.track <- 
        #   track.list[[this.run]][[length(track.list[[this.run]])]]$progress.track <- 
        progress.track
      # name the iterations
      #names(track.list[[this.run]]) <- c("init", paste0("iter", iter))
    }
  }
  names(hclust.list) <- paste0("run", 1:runs)
  # if(verbose){
  #   # names(track.list) <- paste0("run", 1:runs)
  #   out <- list(hclust.list = hclust.list,
  #               track.list = track.list)
  # } else {
  if(runs > 1){
    out <- hclust.list
  } else {
    out <- hclust.list[[1]]
  }
  # }
  invisible(out)
}

#' b-cluster analysis by non-hierarchical iterative ascent clustering
#' strategy
#' 
#' Non-hierarchical b-cluster analysis transfers assessors iteratively to
#' reach a local maximum in sensory differentiation retained.
#' @name bcluster.n
#' @aliases bcluster.n
#' @usage bcluster.n(X, G, M = NULL, measure = "b", max.iter = 500, runs = 1,
#' X.input = "data", tol = exp(-32), seed = 2021)
#' @param X CATA data organized in a three-way array (assessors, products, 
#' attributes)
#' @param G number of clusters (required for non-hierarchical algorithm)
#' @param M initial cluster memberships (default: \code{NULL}), but can be a vector
#' (one run) or a matrix (consumers in rows; runs in columns)
#' @param measure \code{b} (default) for the \code{b}-measure is implemented
#' @param max.iter maximum number of iteration allowed (default \code{500})
#' @param runs number of runs (defaults to \code{1})
#' @param X.input either \code{"data"} (default) or \code{"bc"} if \code{X} is
#' obtained from the function \code{\link[cata]{barray}} 
#' @param tol algorithm stops if variance over 5 iterations is less than 
#' \code{tol} (default: \code{exp(-32)})
#' @param seed for reproducibility (default is \code{2021})
#' @return An object of class \code{bclust.n} (or a list of such objects 
#' if \code{runs>1}), where each such object has the following components: 
#' \itemize{
#' \item{\code{cluster} : vector of the final cluster memberships}
#' \item{\code{totalB} : value of the total sensory differentiation in data set}
#' \item{\code{retainedB} : value of sensory differentiation retained in b-cluster
#' analysis solution}
#' \item{\code{progression} : vector of sensory differentiation retained in each 
#' iteration}
#' \item{\code{iter} : number of iterations completed}
#' \item{\code{finished} : boolean indicates whether the algorithm converged 
#' before \code{max.iter}}}
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Meyners, M., Varela, P., & Næs, T. (2022). 
#' Clustering consumers based on product discrimination in check-all-that-apply 
#' (CATA) data. \emph{Food Quality and Preference}, 104564. 
#' \doi{10.1016/j.foodqual.2022.104564}.
#' @examples
#' data(bread)
#' 
#' # b-cluster analysis on the first 8 consumers and the first 5 attributes
#' (b <- bcluster.n(bread$cata[1:8, , 1:5], G=2))
bcluster.n <- function(X, G, M = NULL, measure = "b", max.iter = 500, runs = 1,
                       X.input = "data", tol = exp(-32), seed = 2021){
  bcluster.call <- match.call()
  
  if(!is.null(M) & runs > 1){
    print("Cluster memberships in M will be the starting point for run 1.")
    print("Subsequent runs will use random start.")
    user.yn <- readline("Continue? (y/n): ")
    if(user.yn != "y"){
      return(print("Cluster analysis stopped"))
    }
  } 
  if(!(measure %in% c("b", "Q"))){
    return(print("Value for measure must be one of \"b\" or \"Q\"", 
                 quote = FALSE))
  }
  
  ## Functions
  getQ <- function(X){
    if(length(dim(X))<3){
      #out <- cochranQ(X, quiet = TRUE)
      return(cochranQ(X, quiet = TRUE))
    } else {
      #out <- sum(apply(X, 3, cochranQ, quiet = TRUE), na.rm=TRUE)
      return(sum(apply(X, 3, cochranQ, quiet = TRUE), na.rm=TRUE))
    }
    #return(out)
  }
  
  # .getCurrent get sensory differentiation retained based on current clusters
  .getCurrent <- function(M, X = NULL, 
                          G = length(unique(M)), 
                          X.bc = NULL, X.bc.dim = NULL, measure = "b"){
    current <- rep(NA, G)
    if(measure == "b"){
      if(!all.equal(dim(X.bc), X.bc.dim)){
        X.bc <- array(X.bc, X.bc.dim)
      }  
    } 
    oneM <- (X.bc.dim[length(X.bc.dim)]==1)
    b.bygroup <- function(g, M, oneM){
      g.mem <- which(M==g)
      return(getb(X.bc[g.mem,,1,], X.bc[g.mem,,2,],
                  oneI = (length(g.mem)==1), oneM = oneM))
    }
    if(measure == "b"){
      return(mapply(b.bygroup, 1:G, MoreArgs = list(M=M, oneM=oneM)))
    }
    if(measure == "Q"){
      for(g in 1:G){
        g.mem <- which(M==g)
        current[g] <- getQ(X[g.mem,,])
      }
    } 
  }
  
  # used in .initCandidate() to evaluate b-measure of cluster g where 'assessor'
  # is either removed (new.g is NA) or added (to cluster number 'new.g')
  .getCandidateMatrix <- function(g, new.g, assessor, M=M, X.bc=X.bc, X.bc.dim=X.bc.dim){
    this.M <- M
    this.M[assessor] <- new.g
    g.mem <- which(this.M==g)
    return(
      getb(X.bc[g.mem,,1,],
           X.bc[g.mem,,2,],
           oneI = (length(g.mem)==1),
           oneM = (X.bc.dim[length(X.bc.dim)]==1))
    )
  }
  
  # .initCandidate calculates the "calculated values matrix"
  .initCandidate <- function(M, X = NULL, G = length(unique(M)), 
                             X.bc = NULL, X.bc.dim = NULL, measure = "b"){
    if(measure %in% "b"){
      # columns: current group, potential group, change in b
      indxGrid <- cbind(curr.g = rep(M, times = G), # 1 
                        change.curr.g = rep(1:G, each = length(M)), # 2
                        new.g = rep(1:G, each = length(M)), # 3
                        assessor = 1:length(M), # 4
                        b = NA) # 5
      # dont' need to calculate if the group doesn't change
      indxGrid[which(indxGrid[,1]==indxGrid[,2]), 3] <- NA
        return(
          matrix(mapply(.getCandidateMatrix, 
                        g=indxGrid[,2], 
                        new.g=indxGrid[,3],
                        assessor=indxGrid[,4],
                        MoreArgs = list(M=M, X.bc=X.bc, X.bc.dim=X.bc.dim)),
                 nrow=length(M), ncol = G))
    } else {
      for(i in seq_along(M)){
        this.M <- M
        curr.g <- M[i]
        for(g in 1:G){
          # drop i from g if in the group, add i to g if not in the group
          this.M[i] <- ifelse(g == curr.g, NA, g)
          if(measure == "Q"){
            out[i,g] <- getQ(X[which(this.M==g),,])
          }
        }
      }
      return(out)
    }
  }
  
  # .updateCandidate updates the "calculated values matrix" 
  .updateCandidate <- function(M, X, candidate, update.i, update.g, 
                               G=length(unique(M)), 
                               X.bc = NULL, X.bc.dim = NULL, measure = "b"){
    if(!all.equal(dim(candidate), c(length(M), G))) return(print("error"))
    if(measure == "b"){
      if(G <= 2){
        return(
          .initCandidate(M, G = length(unique(M)), 
                         X.bc = X.bc, X.bc.dim = NULL, measure = "b")
        )
      } else {
        out <- candidate
        
        # update groups in update.g (all consumers)
        for(g in update.g){
          for(i in seq_along(M)){
            this.M <- M
            curr.g <- M[i]
            this.M[i] <- ifelse(g == curr.g, NA, g) 
            g.mem <- which(this.M==g)
            out[i,g] <- getb(X.bc[g.mem,,1,], 
                             X.bc[g.mem,,2,], 
                             oneI = (length(g.mem)==1),
                             oneM = (X.bc.dim[length(X.bc.dim)]==1)) 
          }
        }
        # update consumers in update.i (groups in update.g already done)
        for(i in update.i){
          for(g in 1:G){
            if(!(g %in% update.g)){
              this.M <- M
              curr.g <- M[i]
              this.M[i] <- ifelse(g == curr.g, NA, g) 
              g.mem <- which(this.M==g)
              out[i,g] <- getb(X.bc[g.mem,,1,], 
                               X.bc[g.mem,,2,], 
                               oneI = (length(g.mem)==1),
                               oneM = (X.bc.dim[length(X.bc.dim)]==1)) 
            }
          }
        }
        # for(i in seq_along(M)){
        #   this.M <- M
        #   curr.g <- M[i]
        #   for(g in 1:G){
        #     if((i %in% update.i) || (g %in% update.g)){
        #       # drop i from g if in the group, add i to g if not in the group
        #       this.M[i] <- ifelse(g == curr.g, NA, g) 
        #       g.mem <- which(this.M==g)
        #       out[i,g] <- getb(X.bc[g.mem,,1,], 
        #                          X.bc[g.mem,,2,], 
        #                          oneI = (length(g.mem)==1),
        #                          oneM = (X.bc.dim[length(X.bc.dim)]==1)) 
        #     }
        #   }
        # }
        return(out)
      }
    } else {
      out <- candidate
      for(i in seq_along(M)){
        this.M <- M
        curr.g <- M[i]
        for(g in 1:G){
          if((i %in% update.i) || (g %in% update.g)){
            # drop i from g if in the group, add i to g if not in the group
            this.M[i] <- ifelse(g == curr.g, NA, g) 
            out[i,g] <- getQ(X[which(this.M==g),,])
          }
        }
      }
      return(out)
    }
  }
  
  .getGainMatrix <- function(g, new.g, assessor, M=M, X.bc=X.bc, X.bc.dim=X.bc.dim){
    this.M <- M
    this.M[assessor] <- new.g
    g.mem <- which(this.M==g)
    return(
      getb(X.bc[g.mem,,1,], X.bc[g.mem,,2,],
           oneI = (length(g.mem)==1),
           oneM = (X.bc.dim[length(X.bc.dim)]==1))
    )
  }
  
  # helper function for .getGainMat
  .thisgain <- function(i, g, current = current, candidate = candidate, M = M){
    if(M[i] == g){
      return(NA)
    } else {
      return(sum(candidate[i, c(g, M[i])]) - sum(current[c(g, M[i])]))
    } 
  }

  # .getGainMat calculates the update (gain) matrix U
  .getGainMat <- function(current, candidate, M, G=length(unique(M))){
      tmp <- cbind(assessor = rep(seq_along(M), times = G), 
                   group = rep(1:G, each = length(M)),
                   gain = NA)
      return(matrix(mapply(.thisgain, 
                           i = rep(seq_along(M), times = G), 
                           g = rep(1:G, each = length(M)), 
                           MoreArgs = list(current = current,
                                           candidate = candidate,
                                           M = M)), nrow = length(M)))
  }
  
  .idSwap <- function(current, candidate, M, G=length(unique(M))){
    gainMat <- .getGainMat(current, candidate, M, G=length(unique(M)))
    best.change <- max(gainMat, na.rm = TRUE)
    if(best.change<0){
      return(list(i=0, g=0))
    } else {
      best.indx <- which(gainMat==best.change, arr.ind = TRUE)
      if(dim(best.indx)[1] == 1){
        return(list(i=best.indx[1], g=best.indx[2]))  
      } else {
        # if there are multiple, pick the best one
        best.indx.row <- sample(dim(best.indx)[1], 1)
        return(list(i=best.indx[best.indx.row, 1], 
                    g=best.indx[best.indx.row, 2]))  
      }
    }
  }
  
  findJ <- function(njj){
    # we know that j is larger than 1 and less than njj
    j.guess <- 1
    njj.sum <- 0
    continue <- TRUE
    while(continue){
      if(njj.sum == njj){
        continue <- FALSE
      } else {
        njj.sum <- njj.sum + j.guess
        j.guess <- j.guess + 1
      }
      if(njj.sum > njj){
        j.guess <- NA
        continue <- FALSE
      }
    }
    return(j.guess)
  }
  
  # get M for a single run
  .getM <- function(nI, G, seed){
    set.seed(seed)
    return(c(sample(G, G, replace = FALSE), sample(G, nI-G, replace = TRUE)))
  }
  
  # get Ms for multiple runs
  startMs <- function(nI, G, Seeds){
    if(!is.null(M)){
      return(cbind(M, mapply(.getM, rep(nI, runs-1), rep(G, runs-1), 
                             (Seeds+1)[-1])))
    } else {
      return(mapply(.getM, rep(nI, runs), rep(G, runs), Seeds+1))
    }
  }
  
  ## End of functions
  
  if(measure %in% "b" && X.input == "data"){
    if(length(dim(X)) == 2){
      X <- array(X, c(dim(X)[1], dim(X)[2], 1), 
                 dimnames = 
                   list(dimnames(X)[[1]], dimnames(X)[[2]], "data"))
    }
    nI <- dim(X)[1]
    nJ <- dim(X)[2]
    nM <- dim(X)[3]
    X.bc <- barray(X)
    nJJ <- dim(X.bc)[2]
  }
  if(measure == "b" && is.null(X.bc)){
    X.bc <- barray(X)
    if(!is.null(X.bc.dim[1])){
      if(!all.equal(dim(X.bc), X.bc.dim)){
        X.bc <- array(X.bc, X.bc.dim)
      }
    }  
  } 
  if (X.input %in% "bc"){
    if(measure %in% "Q"){
      return(print("b-cluster analysis requires raw data for measure Q"))
    }
    X.bc <- X
  }
  msg <- NULL
  if(length(dim(X.bc)) == 3){
    msg <- "Cluster analysis will proceed assuming there is only 1 response variable"
    if(!is.null(msg)) print(msg)
    X.bc <- array(X.bc, c(dim(X.bc)[1], dim(X.bc)[2], dim(X.bc)[3], 1), 
               dimnames = 
                 list(dimnames(X.bc)[[1]], dimnames(X.bc)[[2]], 
                      dimnames(X.bc)[[3]], "data"))
  }
  if (X.input %in% "bc"){
    # X.bc <- X
    X <- NULL
    nI <- dim(X.bc)[1]
    nJJ <- dim(X.bc)[2]
    nM <- dim(X.bc)[4] # [5]
    nJ <- findJ(dim(X.bc)[2])
  }
  
  if(is.null(dimnames(X.bc)[[1]])) dimnames(X.bc)[[1]] <- 1:nI
  if(is.null(dimnames(X.bc)[[2]])) dimnames(X.bc)[[2]] <- 1:nJJ
  if(is.null(dimnames(X.bc)[[3]])) dimnames(X.bc)[[3]] <- letters[2:3]
  if(is.null(dimnames(X.bc)[[4]])) dimnames(X.bc)[[4]] <- 1:nM
  X.bc.dim <- c(nI, nJJ, 2, nM)
  
  if(measure %in% "b"){
    totalB <- sum(abs(X.bc[,,1,] - X.bc[,,2,]))
  }
  if(measure %in% "Q"){
    totalQ <- sum(apply(X, c(1,3), stats::sd) > 0) * (nJ-1)
  }

  .bclustn <- function(M, X.bc, X.bc.dim, measure = "b"){
    current <- .getCurrent(M = M, X.bc = X.bc, X.bc.dim = X.bc.dim, 
                           measure=measure)
    candidate <- .initCandidate(M, X = NULL, X.bc = X.bc, X.bc.dim = X.bc.dim, 
                                measure = measure)
    it <- 0
    continue <- TRUE
    len.Qual <- 1
    Qual <- c(sum(current), rep(NA, max.iter))
    while (continue){
      this.swap <- .idSwap(current, candidate, M)
      if(this.swap$i > 0){
        old.g <- M[this.swap$i]
        # Update group membership vector
        M[this.swap$i] <- this.swap$g
        # Update candidate matrix
        candidate <- .updateCandidate(M, X = NULL, candidate, 
                                      update.i = this.swap$i, 
                                      update.g = c(old.g, this.swap$g), 
                                      G=G, X.bc = X.bc, X.bc.dim = X.bc.dim,
                                      measure = measure)
        # Recalculate current
        current <- .getCurrent(M = M, X.bc = X.bc, 
                               X.bc.dim = X.bc.dim, measure=measure)
        # Update quality measure
        Qual[it+2] <- sum(current) 
        len.Qual <- len.Qual+1
      } else {
        continue <- FALSE # no further swaps improve the solution
        break
      }
      it <- it + 1
      if(it > max.iter){
        continue <- FALSE
        break
      }
      if(len.Qual > 6){
        if (Qual[len.Qual] - Qual[len.Qual-5] < tol){
          continue <- FALSE
          break
        } 
      } 
    }
    o <- 
      list(cluster = M, totalB = sum(abs(X.bc[,,1,]-X.bc[,,2,])),
           retainedB = sum(current), 
           progression = Qual[1:len.Qual], iter = it,
           finish = it<=max.iter)
    class(o) <- "bclust.n"
    return(o)
  }
  
  if(runs > 1){
    if(is.null(M[[1]])){
      Ms <- startMs(nI, G, seed:(seed+runs-1))
    } 
    Ls <- lapply(seq_len(ncol(Ms)), function(i) Ms[,i])
    invisible(lapply(Ls, function(m){
      .bclustn(M = unlist(m), X.bc = X.bc, X.bc.dim = X.bc.dim, measure = measure)
    }))
  } else {
    if(is.list(M)){ 
      M <- unlist(M, use.names = FALSE)
    } else {
      if(is.null(M[[1]])){
        M <- .getM(nI, G, seed)
      }
    }
    invisible(.bclustn(M = unlist(M), X.bc = X.bc, X.bc.dim = X.bc.dim, measure = measure))
  }
}

#' Inspect/summarize many b-cluster analysis runs
#'
#' Inspect many runs of b-cluster analysis. Calculate sensory differentiation 
#' retained and recurrence rate.  
#' @name inspect
#' @aliases inspect
#' @usage inspect(X, G = 2, bestB = NULL, bestM = NULL, inspect.plot = TRUE)
#' @param X list of multiple runs of b-cluster analysis results from 
#' \code{\link[cata]{bcluster.n}} or \code{\link[cata]{bcluster.h}} 
#' @param G number of clusters (required for non-hierarchical algorithm)
#' @param bestB  total sensory differentiation retained in the best solution. If
#' not provided, then \code{bestB} is determined from best solution in the runs
#' provided (in \code{X}).
#' @param bestM  cluster memberships for best solution. If not provided, then 
#' the best solution is determined from the runs provided (in \code{X}).
#' @param inspect.plot default (\code{TRUE}) plots results from the 
#' \code{\link[cata]{inspect}} function 
#' @return A data frame with unique solutions in rows and the following
#' columns:
#' \itemize{
#' \item{\code{B} : Sensory differentiation retained}
#' \item{\code{PctB} : Percentage of the total sensory differentiation retained}
#' \item{\code{B.prop} : Proportion of sensory differentiation retained compared
#' to best solution}
#' \item{\code{Raw.agree} : raw agreement with best solution}
#' \item{\code{Count} : number of runs for which this solution was observed} 
#' \item{\code{Index} : list index (i.e., run number) of first solution  
#' solution in \code{X} corresponding to this row}}
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Meyners, M., Varela, P., & Næs, T. (2022). 
#' Clustering consumers based on product discrimination in check-all-that-apply 
#' (CATA) data. \emph{Food Quality and Preference}, 104564. 
#' \doi{10.1016/j.foodqual.2022.104564}.
#' @examples
#' data(bread)
#' 
#' res <- bcluster.n(bread$cata[1:8, , 1:5], G = 3, runs = 3)
#' (ires <- inspect(res))
#' # get index of solution retaining the most sensory differentiation (in these runs)
#' indx <- ires$Index[1]
#' # cluster memberships for solution of this solution
#' res[[indx]]$cluster
inspect <- function(X, G = 2, bestB = NULL, bestM = NULL, inspect.plot = TRUE){
  # Functions
  get.retainedB <- function(y, k){
    lost.B <- sum(y$height[1:(length(y$height)-k+1)])
    retained.B <- rev(y$retainedB)[k]
    tot.B <- retained.B + lost.B
    return(c(retained.B, tot.B))
  }
  get.raw.agreement <- function(x, y){
    tblx <- table(x)
    tblr <- table(y)
    tblxr <- table(x,y)
    return(1 + (sum(tblxr^2) - sum(tblx^2)/2 - sum(tblr^2)/2) /
             choose(length(x), 2))
  }
  # End functions
  n.sol <- length(X)
  if(!is.null(bestM[1])){
    if(length(bestM) != n.sol){
      return(print("Number of consumers in X must be equal to number of cluster  
                   memberships in bestM"))
    }
  }
  continue <- inspection.complete <- FALSE
  if(class(X[[1]]) %in% "hclust"){
    nI <- nrow(X[[1]]$merge) + 1
    # Mmat   :  membership matrix: rows=assessors, columns=runs
    Mmat <- matrix(unlist(lapply(X, stats::cutree, k=G)), nrow = nI, byrow = FALSE)
    # bb     :  sensory differentiation retained
    bb <- matrix(unlist(lapply(X, get.retainedB, G)), nrow = 2, byrow = FALSE)
    continue <- TRUE
  }
  if (class(X[[1]]) %in% "bclust.n"){
    nI <- length(X[[1]]$cluster)
    # Mmat   :  membership matrix: rows=assessors, columns=runs
    Mmat <- matrix(unlist(lapply(X, function(x){
      return(x$cluster)})), nrow = nI, byrow = FALSE)
    # bb     :  sensory differentiation retained
    bb <- matrix(unlist(lapply(X, function(x){ return(c(x$retainedB,
                                                        x$totalB))})), nrow = 2)
    continue <- TRUE
  }
  if(continue){
    if(is.null(bestB)) bestB <- max(bb[1,])
    maxB <- max(bb)
    # bbMmat :  both, rows=runs, columns=assessors
    bbMmat <- cbind(rep(1, n.sol), bb[1,], t(Mmat))
    # best.sol.indx :  row index of bbMmat having the best solution
    best.sol.indx <- which(bbMmat[,2] == max(bbMmat[,2]))
    # all.sol :  unique solutions
    all.sol <- stats::aggregate(bbMmat[,1] ~ ., data = bbMmat, sum)[,-1]
    # all.sol :  order from best to worst
    all.sol <- all.sol[order(all.sol[,1], decreasing = TRUE), ]
    all.sol.save <- all.sol
    # if B is identical in two or more rows, merge if they agree perfectly
    if(nrow(all.sol) > 1){
      for(rr in nrow(all.sol):2){
        # compare with all rows above with equal B (not just row above)
        rr2.indx <- which(all.sol[1:(rr-1),1] == all.sol[rr,1])
        if(length(rr2.indx)>0){
          for(rr2 in rr2.indx){
            if(isTRUE(all.equal(all.sol[rr,1], all.sol[rr2,1]))){
              #print(paste("rr-1 is the same as rr =", rr))
              # identical B so check if agreement is perfect
              if(isTRUE(all.equal(get.raw.agreement(
                unlist(all.sol[rr, -c(1,ncol(all.sol))]),
                unlist(all.sol[rr2, -c(1,ncol(all.sol))])), 1))){
                #print(paste("raw agreement for rr-1 and rr =", rr))
                all.sol[rr2, ncol(all.sol)] <- 
                  all.sol[rr2, ncol(all.sol)] + all.sol[rr, ncol(all.sol)]
                all.sol[rr, ncol(all.sol)] <- 0
              }
            }
          }
        }
      }
      redundant.indx <- which(all.sol[, ncol(all.sol)] == 0)
      if(length(redundant.indx)>0){
        all.sol <- all.sol[-which(all.sol[, ncol(all.sol)] == 0), ]
      }
    }
    # order by B, then by count
    all.sol <- all.sol[order(all.sol[,1], 
                             all.sol[, ncol(all.sol)], decreasing = TRUE),]
    # best.sol:  cluster memberships of best solution
    best.sol <- unlist(all.sol[1,-c(1,ncol(all.sol))])
    # this step can be improved later
    if(is.null(bestM)) bestM <- best.sol
    raw.a <- apply(as.matrix(all.sol[,-c(1,ncol(all.sol))]), 1,
                   function(x, ref=bestM){
                     tblx <- table(x)
                     tblr <- table(ref)
                     tblxr <- table(x,ref)
                     return(1 + (sum(tblxr^2) -
                                   sum(tblx^2)/2 -
                                   sum(tblr^2)/2) /
                              choose(length(x), 2)) } )
    all.solutions <-
      cbind(all.sol[,1], 
            100*(round(all.sol[,1]/max(bb)[1], 4)),
            all.sol[,1]/bestB, 
            round(raw.a,4), all.sol[,ncol(all.sol)], 
            as.numeric(rownames(all.sol)),
            as.matrix(all.sol[,-c(1,ncol(all.sol))]))
    colnames(all.solutions)[1:6] <- c("B", "PctB", "B.prop", "Raw.agree", 
                                      "Count", "Index")
    colnames(all.solutions)[-c(1:6)] <- paste0("c.", 1:nI)
    rownames(all.solutions) <- NULL
    if(inspect.plot){
      curr.par <- graphics::par(no.readonly = TRUE) # current parameters
      on.exit(graphics::par(curr.par)) # return to current parameters on exit
      graphics::par(mar=c(7, 4, 4, 2) + 0.1)
      plot(c(1,0), range(all.sol[,1]/bestB), type="n", xlim = c(1,0),
           main = paste0("Strawberry cluster analysis (G=", G, ") from ", 
                         n.sol, " runs"),
           ylab = "Sensory Differentiation Retained (B)",
           xlab = "Proportion of runs with this B or higher",
           sub = paste0("Labels indicate % agreement with best solution"))
      graphics::text(x = cumsum(all.sol[,ncol(all.sol)]) /
                       sum(all.sol[,ncol(all.sol)]),
                     y = all.sol[,1]/bestB,
                     labels = as.character(round(raw.a,2)*100), srt=45)
    }
    inspection.complete <- TRUE
  }
  if(!inspection.complete){
    return("Input must be a list of b-cluster analysis results (many runs)")
  } else {
    num.best <- length(all.solutions[,1] == all.solutions[1,1])
    if(num.best > 1){
      cat("The best solution from these ", n.sol, " runs (row 1) ",
      "is compared with the solutions found in other runs.", 
      "\n", sep = "")
    }
    all.solutions <- as.data.frame(all.solutions)
    all.solutions[,6] <- as.integer(all.solutions[,6])
    return(all.solutions[,1:6])
  }
} 

#' Calculate the b-measure
#'
#' Function to calculate the b-measure, which quantifies the sensory 
#' differentiation retained. 
#' @name getb
#' @aliases getb
#' @usage getb(X.b, X.c, oneI = FALSE, oneM = FALSE)
#' @param X.b three-way (\code{I, J(J-1)/2, M}) array with 
#' \code{I} assessors, \code{J(J-1)/2} product comparisons, \code{M} CATA 
#' attributes, where values are counts of type \code{b} from the function 
#' \code{\link[cata]{barray}})
#' @param X.c array of same dimension as \code{X.b}, where values are counts of 
#' type \code{b} from the function \code{\link[cata]{barray}})
#' @param oneI indicates whether calculation is for one assessor (default: 
#' \code{FALSE})
#' @param oneM  indicates whether calculation is for one attribute (default: 
#' \code{FALSE})
#' @return b-measure
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Meyners, M., Varela, P., & Næs, T. (2022). 
#' Clustering consumers based on product discrimination in check-all-that-apply 
#' (CATA) data. \emph{Food Quality and Preference}, 104564. 
#' \doi{10.1016/j.foodqual.2022.104564}.
#' @examples
#' data(bread)
#' 
#' bread.bc <- barray(bread$cata[1:8,,1:5])
#' getb(bread.bc[,,1,], bread.bc[,,2,])
getb <- function(X.b, X.c, oneI = FALSE, oneM = FALSE){
  if(any(oneI, oneM)) {
    if(all(oneI, oneM)){
      dims <- c(1, length(X.b), 1)
    } else {
      dims <- dim(X.b) # 2D
      if(oneI){
        dims <- c(1, dims) # nI was 1
      }
      if(oneM){
        dims <- c(dims, 1) # nM was 1
      }
    }
    # make 3-way array
    X.b <- array(X.b, dims)
    X.c <- array(X.c, dims)
  }
  dims <- dim(X.b) # dimension is now nI x nJJ x nM
  if(!all.equal(dims, dim(X.c))){
    # dimension should also be nI x nJJ x nM
    return("Cannot calculate b-measure for unequal sized arrays")
  }
  # keep 3-way array structure
  return(sum(((apply(array(X.b - X.c, dims), 2:3, sum)^2) /
                (apply(array(X.b + X.c, dims), 2:3, sum))), na.rm = TRUE))
  
}

#' Cochran's Q test
#'
#' Calculate Cochran's Q test statistic. The null hypothesis that is assumed is 
#' that product proportions are all equal. The alternative hypothesis is that 
#' product proportions are not all equal.
#' @name cochranQ
#' @aliases cochranQ
#' @usage cochranQ(X, na.rm = TRUE, quiet = FALSE, digits = getOption("digits"))
#' @param X matrix of \code{I} assessors (rows) and \code{J} products (columns)
#' where values are \code{0} (not checked) or \code{1} (checked)
#' @param na.rm should \code{NA} values be removed?
#' @param quiet if \code{FALSE} (default) then it prints information related to 
#' the test; if \code{TRUE} it returns only the test statistic (\code{Q})
#' @param digits significant digits (to display)
#' @return Q test statistic
#' @export
#' @encoding UTF-8
#' @seealso \code{\link[cata]{mcnemarQ}}
#' @references  
#' 
#' Cochran, W. G. (1950). The comparison of percentages in matched samples. 
#' \emph{Biometrika}, 37, 256-266.
#' 
#' Meyners, M., Castura, J.C., & Carr, B.T. (2013). Existing and  
#' new approaches for the analysis of CATA data. \emph{Food Quality and Preference}, 
#' 30, 309-319, \doi{10.1016/j.foodqual.2013.06.010}
#' @examples
#' data(bread)
#' 
#' # Cochran's Q test on the first 25 consumers on the first attribute ("Fresh")
#' cochranQ(bread$cata[1:25,,1])
cochranQ <- function(X, na.rm = TRUE, quiet = FALSE, 
                     digits = getOption("digits")){
  if(is.vector(X)){
    X <- matrix(X, nrow = 1)
  }
  if(na.rm){
    X <- X[stats::complete.cases(X),]
    if(is.vector(X)){
      X <- matrix(X, nrow = 1)
    }
  }
  X <- X[which(rowSums(X) < ncol(X)),]
  if(is.vector(X)){
    X <- matrix(X, nrow = 1)
  }
  if(sum(X)==0){
    return(Q = 0)
  } 
  numQ <- (ncol(X) * (ncol(X)-1) * (sum(colSums(X)^2) - 
                                      (sum(colSums(X))^2)/ncol(X)))
  denomQ <-  ncol(X) * sum(rowSums(X))-sum(rowSums(X)^2)
  Q = numQ / denomQ
  if(!quiet){
    cat(paste(c("",
                "Cochran's Q test", 
                "----------------", "",
                "H0: Citation rates equal for all products (columns)",
                "H1: Citation rates are not equal for all products (columns)",
                ""), 
              collapse = '\n'))
    print(data.frame(Q  = signif(Q, digits = digits), 
                     df = ncol(X)-1,
                     p.value = 
                       signif(stats::pchisq(Q, ncol(X)-1, lower.tail = FALSE),
                              digits = digits),
                     row.names = ""))
  }
  invisible(Q)
}

#' McNemar's test
#'
#' Pairwise tests are conducted using the two-tailed binomial test. These tests
#' can be conducted after Cochran's Q test. 
#' @name mcnemarQ
#' @aliases mcnemarQ
#' @usage mcnemarQ(X, na.rm = TRUE, quiet = FALSE, digits = getOption("digits"))
#' @param X matrix of \code{I} assessors (rows) and \code{J} products (columns)
#' where values are \code{0} (not checked) or \code{1} (checked)
#' @param na.rm should \code{NA} values be removed?
#' @param quiet if \code{FALSE} (default) then it prints information related to 
#' the test; if \code{TRUE} it returns only the test statistic (\code{Q})
#' @param digits significant digits (to display)
#' @return Test results for all McNemar pairwise tests conducted via the 
#' binomial test
#' @export
#' @encoding UTF-8
#' @seealso \code{\link[cata]{cochranQ}}
#' @references  
#' 
#' Cochran, W. G. (1950). The comparison of percentages in matched samples. 
#' \emph{Biometrika}, 37, 256-266. 
#' 
#' McNemar, Q. (1947). Note on the sampling error of the difference between 
#' correlated proportions or percentages. \emph{Psychometrika}, 12(2), 153-157.
#' 
#' Meyners, M., Castura, J.C., & Carr, B.T. (2013). Existing and  
#' new approaches for the analysis of CATA data. \emph{Food Quality and 
#' Preference}, 30, 309-319, \doi{10.1016/j.foodqual.2013.06.010}
#' 
#' @examples
#' data(bread)
#' 
#' # McNemar's exact pairwise test for all product pairs
#' # on the first 25 consumers and the first attribute ("Fresh")
#' mcnemarQ(bread$cata[1:25,,1])
mcnemarQ <- function(X, na.rm = TRUE, quiet = FALSE, 
                     digits = getOption("digits")){
  if(is.vector(X)){
    return(print(
      "Input X must be a matrix with at least 2 rows and 2 columns"))
  }
  if(na.rm){
    X <- X[stats::complete.cases(X),]
    if(is.vector(X)){
      return(print(
        "Input X must be a matrix with at least 2 rows and 2 columns"))
    }
  }
  X <- X[which(rowSums(X) < ncol(X)),]
  if(is.vector(X)){
    X <- matrix(X, nrow = 1)
  }
  if(any(dim(X) == 1)){
    return(print(
      "Input X must be a matrix with at least 2 rows and 2 columns"))
  }
  if(sum(X)==0 || sum(X)==prod(dim(X))){
    return(print("Input X must be a matrix with variability"))
  }
  res <- matrix(NA, nrow = choose(ncol(X), 2), ncol = 6,
                dimnames = list(NULL, c("index.1", "index.2", "b", "c", 
                                        "proportion", "p.value")))
  res[,1:2] <- t(utils::combn(ncol(X), 2))
  res[,3:4] <- apply(barray(X), 2:3, sum)
  res[,5] <- res[,3] / rowSums(res[,3:4])
  res[,6] <- mapply(function(x, y){ 
    2*sum(stats::dbinom(0:min(x, sum(x,y)-x), size = sum(x,y), 
                        prob = .5))}, res[,3], res[,4])
  res[,6] <- mapply(min, res[,6], 1)

  if(!quiet){
    cat(paste(c("",
                "McNemar's pairwise test", 
                "-----------------------", "",
                "H0: Citation rates equal for both products",
                "H1: Citation rates are not equal for both products",
                "", ""), 
              collapse = '\n'))
    res.print <- res
    res.print[,5:6] <- signif(res.print[,5:6], digits = digits)
    print(as.data.frame(res.print, row.names = ""))
  }
  invisible(res)
}

#' Convert 3d array of CATA data to 4d array of CATA differences
#'
#' Converts a three-dimensional array (\code{I} assessors, \code{J}
#' products, \code{M} attributes) to a four-dimensional array of product
#' comparisons (\code{I} assessors, \code{J(J-1)/2}
#' product comparisons, two outcomes (of type \code{b} or \code{c}), \code{M} 
#' attributes)
#'  
#' @name barray
#' @aliases barray
#' @usage barray(X, values = "bc", type.in = "binary", type.out = "binary")
#' @param X three-dimensional array (\code{I} assessors, \code{J}
#' products, \code{M} attributes) where values are \code{0} (not checked) 
#' or \code{1} (checked)
#' @param values \code{"bc"} (default) returns two outcomes: \code{b} and 
#' \code{c}; otherwise \code{"abcd"} returns four outcomes: \code{a}, \code{b}, 
#' \code{c}, \code{d}.
#' @param type.in type of data submitted; default (\code{binary}) may be set to
#' \code{ordinal} or \code{scale}.
#' @param type.out currently only \code{binary} is implemented
#' @return A four-dimensional array of product comparisons having \code{I} 
#' assessors, \code{J(J-1)/2} product comparisons, outcomes (see \code{values}
#' parameter), \code{M} attributes
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Meyners, M., Varela, P., & Næs, T. (2022). 
#' Clustering consumers based on product discrimination in check-all-that-apply 
#' (CATA) data. \emph{Food Quality and Preference}, 104564. 
#' \doi{10.1016/j.foodqual.2022.104564}.
#' 
#' @examples
#' data(bread)
#' 
#' # Get the 4d array of CATA differences for the first 8 consumers
#' b <- barray(bread$cata[1:8,,])
barray <- function(X, values = "bc", type.in = "binary", type.out = "binary"){
  abcd_1attribute <- function(X, type.in = "binary"){
    .abcd_1pair <- function(x, y, type.in = "binary"){
      x <- c(x) # product 'x' (1 response per consumer)
      y <- c(y) # product 'y' (1 response per consumer)
      if(type.in == "binary"){
        return(list(a = sum(x+y == 2), b = sum(x-y == 1),
                    c = sum(x-y == -1), d = sum(x+y == 0)))
      }
      if(type.in %in% c("scale", "ordinal")){
        return(list(a = sum(all(x == y, x != 0)), 
                    b = sum(x > y),
                    c = sum(x < y), 
                    d = sum(all(x == y, x == 0))))
      }
    }
    if(is.vector(X)){
      X <- matrix(X, nrow = 1)
    } else {
      X <- as.matrix(X) # assessors x products
    }
    pair.items <- utils::combn(ncol(X),2)
    this.pair <- matrix(0, nrow = ncol(pair.items), ncol = 4, 
                        dimnames = list(NULL, letters[1:4]))
    for(i in 1:ncol(pair.items)){
      this.pair[i, ] <- unlist(.abcd_1pair(X[,pair.items[1,i]], 
                                           X[,pair.items[2,i]], type.in = type.in))
    }
    return(this.pair)
  }
  nI <- dim(X)[1]
  nJ <- dim(X)[2]
  if(length(dim(X))==2){
    oX <- X
    X <- array(NA, c(nI, nJ, 1), dimnames = list(dimnames(X)[[1]],
                                                 dimnames(X)[[2]], 1))
    X[,,1] <- oX
  }
  nM <- dim(X)[3]
  out <- array(0, c(nI, nJ*(nJ-1)/2, 4, nM), 
               dimnames = list(dimnames(X)[[1]], 
                               apply(t(utils::combn(nJ,2)), 1, paste0, collapse = "_"), 
                               letters[1:4], 
                               colnames(X[1,,])))
  for (i in 1:nI){
    for(m in 1:nM){
      this.data <- X[i,,m]
      out[i,,1:4,m] <- abcd_1attribute(this.data, type.in = type.in)
    }
  }
  if (values == "bc") {
    return(out[,, 2:3,])
  } else {
    if (values == "abcd") {
      return(out)
    } else {
      return(print(paste("Invalid option: values =", values), quote = FALSE))
    }
  }
}

#' Converts 3d array of CATA data to a wide 2d matrix format
#'
#' Converts a three-dimensional array (\code{I} assessors, \code{J}
#' products, \code{M} attributes) to a two-dimensional matrix
#' (\code{J} products, (\code{I} assessors, \code{M} attributes))
#'  
#' @name toWideMatrix
#' @aliases toWideMatrix
#' @usage toWideMatrix(X)
#' @param X three-dimensional array (\code{I} assessors, \code{J}
#' products, \code{M} attributes) where values are \code{0} (not checked) 
#' or \code{1} (checked)
#' @return A matrix with \code{J} products in rows and \code{I} 
#' assessors * \code{M} attributes in columns
#' @export
#' @encoding UTF-8
#' @examples
#' data(bread)
#' 
#' # convert CATA results from the first 8 consumers and the first 4 attributes
#' # to a wide matrix
#' toWideMatrix(bread$cata[1:8,,1:4])
toWideMatrix <- function(X){
  out <- X[1,,]
  for(i in 2:dim(X)[[1L]]){
    out <- cbind(out, X[i,,])
  }
  return(out)
}

#' Converts 3d array of CATA data to a tall 2d matrix format
#'
#' Converts a three-dimensional array (\code{I} assessors, \code{J}
#' products, \code{M} attributes) to a two-dimensional matrix with
#' (\code{I} assessors, \code{J} products) rows and (\code{M} 
#' attributes) columns, optionally preceded by two columns of row headers.
#'  
#' @name toMatrix
#' @aliases toMatrix
#' @usage toMatrix(X, header.rows = TRUE, oneI = FALSE, oneM = FALSE)
#' @param X three-dimensional array (\code{I} assessors, \code{J}
#' products, \code{M} attributes) where values are \code{0} (not checked) 
#' or \code{1} (checked)
#' @param header.rows \code{TRUE} (default) includes row headers; set to
#' \code{FALSE} to exclude these headers
#' @param oneI indicates whether calculation is for one assessor (default: 
#' \code{FALSE})
#' @param oneM  indicates whether calculation is for one attribute (default: 
#' \code{FALSE})
#' @return A matrix with \code{I} assessors * \code{J} products in rows
#' and \code{M} attributes in columns (preceded by 2 columns)
#' of headers if \code{header.rows = TRUE}
#' @export
#' @encoding UTF-8
#' @examples
#' data(bread)
#' 
#' # convert CATA results from the first 8 consumers and the first 4 attributes
#' # to a tall matrix
#' toMatrix(bread$cata[1:8,,1:4])
toMatrix <- function(X, header.rows = TRUE, oneI = FALSE, oneM = FALSE){
  if(any(oneI, oneM)) {
    if(all(oneI, oneM)){
      dims <- c(1, length(X), 1)
    } else {
      dims <- dim(X) # 2D
      if(oneI){
        dims <- c(1, dims) # nI was 1
      }
      if(oneM){
        dims <- c(dims, 1) # nM was 1
      }
    }
    X <- array(X, c(dims[1], dims[2], dims[3]))
  }
  if(length(dim(X)) == 2){
    X <- array(X, c(dim(X)[1], dim(X)[2], 1), 
               dimnames = list(dimnames(X)[[1]], dimnames(X)[[2]], "response"))
  }
  dims <- dim(X)
  if(length(dims) != 3){
    return("Function is to convert an array to a matrix")
  } 
  nI <- dims[1]
  nJ <- dims[2]
  nM <- dims[3]
  if(is.null(dimnames(X)[[1]])[1]) dimnames(X)[[1]] <- 1:nI
  if(is.null(dimnames(X)[[2]])[1]) dimnames(X)[[2]] <- 1:nJ
  if(is.null(dimnames(X)[[3]])[1]) dimnames(X)[[3]] <- 1:nM
  
  this.rownames <- paste0(rep(dimnames(X)[[1]], each = nJ), "_",
                          rep(dimnames(X)[[2]], times = nI))
  Xout <- matrix(NA, ncol = dim(X)[[3]], nrow = prod(dim(X)[-3]),
                 dimnames = list(this.rownames, dimnames(X)[[3]]))
  for(cc in 1:nM){
    Xout[,cc] <- c(aperm(X[,,cc], 2:1))
  }
  if(header.rows){
    Xout <- data.frame(subject = as.factor(rep(dimnames(X)[[1]], each = nJ)),
                       product = as.factor(rep(dimnames(X)[[2]], times = nI)),
                       Xout)
  }
  return(Xout)
}

#' Calculate RV Coefficient
#'
#' Calculate RV coefficient
#' @name rv.coef
#' @aliases rv.coef
#' @usage rv.coef(X, Y, method = 1)
#' @param X input matrix (same dimensions as \code{Y})
#' @param Y input matrix (same dimensions as \code{X})
#' @param method \code{1} (default) and \code{2} give identical RV coefficients
#' @return RV coefficient
#' @export
#' @encoding UTF-8
#' @references Robert, P., & Escoufier, Y. (1976). A unifying tool for linear 
#' multivariate statistical methods: the RV-coefficient. \emph{Journal of the Royal 
#' Statistical Society: Series C (Applied Statistics)}, 25, 257-265.
#' @examples
#' # Generate some data
#' set.seed(123)
#' X <- matrix(rnorm(8), nrow = 4)
#' Y <- matrix(rnorm(8), nrow = 4)
#' 
#' # get the RV coefficient
#' rv.coef(X, Y)
rv.coef <- function(X, Y, method = 1){
  tr <- function(X){
    sum(diag(X))
  }
  X <- sweep(X, 2, apply(X, 2, mean), "-")
  Y <- sweep(Y, 2, apply(Y, 2, mean), "-")
  XX <- tcrossprod(X,X)
  YY <- tcrossprod(Y,Y)
  if(method==1){
    out <- tr(XX %*% YY) /
      sqrt(tr(XX %*% XX) %*% tr(YY %*% YY) )
  }
  if(method==2){
    out <- sum(c(XX)*c(YY)) /
      sqrt(sum(XX^2)*sum(YY^2))
  }
  return(c(out))
}

#' Salton's cosine measure
#'
#' Calculate Salton's cosine measure
#' @name salton
#' @aliases salton
#' @usage salton(X, Y)
#' @param X input matrix (same dimensions as \code{Y})
#' @param Y input matrix (same dimensions as \code{X})
#' @return Salton's cosine measure
#' @export
#' @encoding UTF-8
#' @references Salton, G., & McGill, M.J. (1983). \emph{Introduction to Modern 
#' Information Retrieval}. Toronto: McGraw-Hill.
#' @examples
#' # Generate some data
#' set.seed(123)
#' X <- matrix(rnorm(8), nrow = 4)
#' Y <- matrix(rnorm(8), nrow = 4)
#' 
#' # get Salton's cosine measure
#' salton(X, Y)
salton <- function(X, Y){
  out <- sum(c(X)*c(Y)) / sqrt(sum(X^2)*sum(Y^2))
  return(out)
}


#' Apply top-k box coding to scale data
#'
#' Apply top-k box coding to scale data. Using defaults give top-2 box (T2B) coding.
#' @name code.topk
#' @aliases code.topk
#' @usage code.topk(X, zero.below = 8, one.above = 7)
#' @param X input matrix
#' @param zero.below default is \code{8}; values below this numeric threshold will be coded \code{0}; use \code{NULL} if there is no such threshold
#' @param one.above default is \code{7}; values above this numeric threshold will be coded \code{1}; use \code{NULL} if there is no such threshold
#' @return matrix \code{X} with top-k coding applied
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Meyners, M., Pohjanheimo, T., Varela, P., & Næs, T. (2023). 
#' An approach for clustering consumers by their top-box and top-choice responses. 
#' \emph{Journal of Sensory Studies}, e12860. \doi{10.1111/joss.12860}
#' @examples
#' # Generate some data
#' set.seed(123)
#' X <- matrix(sample(1:9, 100, replace = TRUE), nrow = 5)
#' 
#' # apply top-2 box (T2B) coding
#' code.topk(X, zero.below = 8, one.above = 7)
code.topk <- function(X, zero.below = 8, one.above = 7){
  Y <- X
  if(!is.null(zero.below) && !is.null(one.above)){
    if(zero.below != one.above + 1){
      return(
        print("zero.below should be equal to one.above+1 if both specified"))
    }
  }
  if(!is.null(zero.below)){
    Y[Y < zero.below] <- 0
  }
  if(!is.null(one.above)){
    Y[Y > one.above] <- 1
  }
  return(Y)
}

#' Apply top-c choices coding to a vector of scale data from a respondent
#'
#' Apply top-c choices coding to a vector of scale data from a respondent
#' @name topc
#' @aliases topc
#' @usage topc(x, c = 2, coding = "B")
#' @param x input matrix
#' @param c number of top choices considered to be 'success'; other choices are 
#' considered to be 'failure' and are coded \code{0}
#' @param coding \code{"B"} (default) codes all successes as \code{1}; 
#' \code{"N"} codes all successes with their numeric coding
#' @return matrix \code{X} with top-k coding applied
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Meyners, M., Pohjanheimo, T., Varela, P., & Næs, T. (2023). 
#' An approach for clustering consumers by their top-box and top-choice responses. 
#' \emph{Journal of Sensory Studies}, e12860. \doi{10.1111/joss.12860}
#' @examples
#' # Generate some data
#' set.seed(123)
#' X <- matrix(sample(1:9, 100, replace = TRUE), nrow = 5)
#' 
#' # apply top-2 choice (T2C) coding
#' apply(X, 1, topc)
topc <- function(x, c = 2, coding = "B"){
  #anything with a value of 0 is special, so give in the max value
  y <- max(x)-x+1
  ranky <- rank(y, ties.method = "average")
  ranky.u <- sort(unique(ranky))
  it <- 0
  indx <- c()
  nc <- 0
  while(nc < c){
    it <- it + 1
    indx <- c(indx, which(ranky == ranky.u[it]))
    nc <- length(indx)
  }
  out <- x*0
  if(coding %in% "B"){
    out[indx] <- 1
  }
  if(coding %in% "N"){
    out[indx] <- x[indx]
  }
  return(out)
}


