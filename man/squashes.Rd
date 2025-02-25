\alias{squashes}
\name{Consumer CATA data set: squashes}
\docType{data}
\title{Consumer CATA data set: squashes}
\description{Check-all-that-apply (CATA) data set from 100 assessors (consumers)
who evaluated 11 products (blackcurrant squashes) using 34 terms.}
\format{
  Array of data coded \code{1} if the assessor checked the term for the attribute
  and \code{0} otherwise.
    \itemize{
      \item{Dimension 1 : 100 assessors}
      \item{Dimension 2 :  11  products}
      \item{Dimension 3 :  34 terms}
    }
  A response is coded \code{1} if the attribute is checked; otherwise it is coded \code{0}.
}
\encoding{UTF-8}
\references{Ng, M., Chaya, C., & Hort, J. (2013). Beyond liking: Comparing the 
measurement of emotional response using EsSense Profile and consumer defined 
check-all-that-apply methodologies. \emph{Food Quality and Preference}, 28, 
193-205, \doi{10.1016/j.foodqual.2012.08.012}}
\examples{
  data(squashes)
  # view CATA frequency table
  apply(squashes, 2:3, sum)
}

