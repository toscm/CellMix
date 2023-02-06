\name{gedAlgorithm.cs_lsfit}
\alias{cs-lsfit-ged}
\alias{gedAlgorithm.cs_lsfit}
\title{Cell-Specific Expression by Standard Least-Squares}
\usage{
  gedAlgorithm.cs_lsfit(..., rescale = TRUE,
    fit = c("ls", "nnls"))
}
\arguments{
  \item{...}{extra arguments passed to fitting the methods
  \code{\link{.nn_lsfit}} or \code{\link[NMF]{.fcnnls}}.}

  \item{rescale}{logical used when estimating proportions
  from signatures, that indicates if the esti,ated
  coefficients should be scaled to sum up to one
  (\code{TRUE}) or left as estimated by the linear
  regression (\code{FALSE}). This scaling is performed
  after the coefficients have been forced to be
  nonnegative.}

  \item{fit}{least-square fitting method: \code{ls} uses
  \code{\link{lm}}, \code{nnls} uses
  \code{\link[NMF]{fcnnls}}.}
}
\description{
  Estimates cell-specific proportions given known
  proportions expression signatures, using least-squares
  fitting.
}
\details{
  The algorithm applies the same methods as the ged
  algorithm \code{\link[=lsfit-ged]{lsfit}} but to the
  transposed problem of estimating signatures from
  proportions. It is included in the \pkg{CellMix} package
  for test/experimental purposes.
}
\examples{
# random target matrix
x <- rmatrix(100, 20)
# random cell proprtions
p <- rmatrix(3, 20)

# deconvolve using standard least-squares
res <- ged(x, p, 'cs-lsfit')
head(basis(res))
# used proportions are stored in the result coefficient matrix
identical(coef(res), p)

# deconvolve using nonnegative least-squares
res <- ged(x, p, 'cs-lsfit', fit = 'nnls')
head(basis(res))
# used proportions are stored in the result coefficient matrix
identical(coef(res), p)
}

