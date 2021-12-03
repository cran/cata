\alias{bread}
\name{Consumer CATA data set: bread}
\docType{data}
\title{Consumer CATA data set: bread}
\description{Raw results from CATA and Liking evaluations of six bread products samples by 161 consumers.}
\format{
A list with 4 items: 
\itemize{
  \item{\code{$cata} :  check-all-that-apply (CATA) data (array, 161 consumers x 6 breads x 31 sensory attributes)}
  \item{\code{$liking} :  9-point hedonic scale data (matrix, 161 consumers x 6 breads)}
  \item{\code{$ideal.cata} :  check-all-that-apply (CATA) data for ideal bread (matrix, 161 consumers x 31 sensory attributes)}
  \item{\code{$liking} :  9-point hedonic scale data for ideal bread(vector, 161 consumers)}
  }

CATA data is coded \code{1} if the attribute is checked; otherwise it is coded \code{0}
  }
\encoding{UTF-8}
\references{Meyners, M., Castura, J.C., & Carr, B.T. (2013). Existing and new approaches for the analysis of CATA data. \emph{Food Quality and Preference}, 30, 309-319, \doi{10.1016/j.foodqual.2013.06.010}}
\examples{
data(bread)
head(bread$cata)
}

