% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cata.R
\encoding{UTF-8}
\name{rv.coef}
\alias{rv.coef}
\title{Calculate RV Coefficient}
\usage{
rv.coef(X, Y, method = 1)
}
\arguments{
\item{X}{input matrix (same dimensions as \code{Y})}

\item{Y}{input matrix (same dimensions as \code{X})}

\item{method}{\code{1} (default) and \code{2} give identical RV coefficients}
}
\value{
RV coefficient
}
\description{
Calculate RV coefficient
}
\examples{
# Generate some data
set.seed(123)
X <- matrix(rnorm(8), nrow = 4)
Y <- matrix(rnorm(8), nrow = 4)

# get the RV coefficient
rv.coef(X, Y)
}
\references{
Robert, P., & Escoufier, Y. (1976). A unifying tool for linear 
multivariate statistical methods: the RV-coefficient. \emph{Journal of the Royal 
Statistical Society: Series C (Applied Statistics)}, 25, 257-265.
}
