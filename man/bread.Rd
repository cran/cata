\alias{bread}
\name{Consumer CATA data set: bread}
\docType{data}
\title{Consumer CATA data set: bread}
\description{Consumer data from 161 assessors (consumers) who evaluated 6 products 
(breads) using a check-all-that-apply (CATA) question and the 9-point hedonic scale.}
\format{
List with items: 
\itemize{
  \item{\code{$cata}:       CATA data array (161 consumers \eqn{\times} 6 breads 
  \eqn{\times} 31 sensory attributes)}
  \item{\code{$liking}:     9-point hedonic scale data matrix (161 consumers \eqn{\times} 6 breads)}
  \item{\code{$ideal.cata}: CATA data matrix for ideal bread evaluations (161 consumers 
  \eqn{\times} 31 sensory attributes)}
  \item{\code{$liking}:     9-point hedonic scale data vector for ideal bread 
  (length: 161 consumers)}
  }
  Responses to CATA questions are each coded \code{1} if the assessor checks the
  attribute for the product; otherwise, it is coded \code{0}
  }
\encoding{UTF-8}
\references{Meyners, M., Castura, J.C., & Carr, B.T. (2013). Existing and new 
approaches for the analysis of CATA data. \emph{Food Quality and Preference}, 
30, 309-319, \doi{10.1016/j.foodqual.2013.06.010}}
\examples{
data(bread)
  # view CATA frequency table
  apply(bread$cata, 2:3, sum)
}

