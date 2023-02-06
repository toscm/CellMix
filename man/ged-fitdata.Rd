\name{basisfit}
\alias{basisfit}
\alias{coeffit}
\title{Accessing Fit Data in Deconvolution Results}
\usage{
  basisfit(x, ...)

  coeffit(x, ...)
}
\arguments{
  \item{x}{an object returned by \code{\link{ged}}, i.e.
  typically an \code{\linkS4class{NMF}} object.}

  \item{...}{extra arguments to allow extension.}
}
\description{
  The function \code{basisfit} is a generic function to
  extract fit data computed when estimating a model's set
  of basis components.

  In the case of gene expression deconvlution, basis
  components consist in the cell-specific -- differential
  -- gene expression values for each estimated cell type.

  The function \code{coeffit} is a generic function to
  extract fit data computed when estimating a model's set
  of coefficients.

  In gene expression deconvolution, these correspond to
  cell type proportions in each sample.
}
\details{
  The default methods for \code{basisfit} (resp.
  \code{coeffit}) returns the "elements" called
  \code{'basisfit'} (resp. \code{'coeffit'}), by as
  returned by \code{x$basisfit} (resp. \code{x$basisfit}).
  This means that they should work on objects for which a
  suitable method \code{$} exists.
}

