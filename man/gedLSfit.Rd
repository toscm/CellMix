\name{.gedLSfit}
\alias{.gedLSfit}
\alias{.nn_lsfit}
\title{Partial Gene Expression Deconvolution by Least-Square}
\usage{
  .gedLSfit(X, seed, rescale = TRUE, fit = c("ls", "nnls"),
    ...)

  .nn_lsfit(x, y, nneg = c("iterate", "pmax", "none"), ...)
}
\arguments{
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

  \item{...}{extra arguments passed to fitting the methods
  \code{\link{.nn_lsfit}} or \code{\link[NMF]{.fcnnls}}.}

  \item{x}{matrix of known cell-specific profiles (i.e.
  cell signatures), with features in rows and cell type in
  columns.}

  \item{y}{matrix of observed mixed expression data, with
  features in rows and samples in columns. The number of
  samples must be greater -- or equal -- than the number of
  cell types.}

  \item{nneg}{specification of the method used to enforce
  the nonnegativity of the estimated proportions.

  Accepted values are: \describe{ \item{\code{'iterate'}:}{
  applies the procedure described in \cite{Abbas2009}. For
  each sample separately, a sequence of least-square fits
  are performed, starting with all cell types, and where
  the cell type corresponding to the lowest negative fitted
  coefficient is excluded from the next fit, and its
  associated final proportion set to zero. This iterative
  process stops when all coefficients are nonnegative.}
  \item{\code{'pmax'}:}{ single least-square fit, where all
  negative estimated proportions are set to zero.}
  \item{\code{NA} or \code{'none'}: }{ single least-square
  fit where the estimated proportions are returned
  unconstrained.} }}

  \item{X}{target matrix}

  \item{seed}{initial NMF model as an \code{\link{NMF}}
  object.}
}
\value{
  an \code{\link{NMF}} object.
}
\description{
  \code{.gedLSfit} provides access to partial deconvolution
  methods that are based on least-squares fits.

  \code{.nn_lsfit} implements a standard least-square fit
  with various different procedures to enforce nonnegative
  coefficients. In particular, it implements the iterative
  procedure described in \cite{Abbas et al. (2009)}.
}
\references{
  Abbas AR, Wolslegel K, Seshasayee D, Modrusan Z and Clark
  HF (2009). "Deconvolution of blood microarray data
  identifies cellular activation patterns in systemic lupus
  erythematosus." _PloS one_, *4*(7), pp. e6098. ISSN
  1932-6203, <URL:
  http://dx.doi.org/10.1371/journal.pone.0006098>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/19568420>.
}
\keyword{internal}

