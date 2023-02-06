\name{gedAlgorithm.lsfit}
\alias{gedAlgorithm.lsfit}
\alias{lsfit-ged}
\title{Partial Gene Expression Deconvolution by Standard Least-Squares}
\usage{
  gedAlgorithm.lsfit(..., rescale = TRUE,
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
  Estimates cell/tissue proportions given a known set of
  cell/tissue-specific expression signatures, using
  standard least-squares as proposed by \cite{Abbas et al.
  (2009)}.
}
\details{
  The default algorithm uses a heuristic to enforce the
  nonnegativity of the estimated proportions, that consists
  in fitting successive regressions, each time excluding
  the most negative coefficient from the model, until all
  coefficients are nonnegative. In this case all
  regressions are fitted using the function
  \code{\link{lm}}.

  An alternative least-square fitting method is included
  for test/experimental purposes. It uses the fast
  combinatorial nonnegative least-square method of
  \cite{Van Benthem et al. (2004)}, which was adapted by
  \cite{Kim et al. (2007)} to perform nonnegative matrix
  factorization of gene expression -- but not originally
  for deconvolution. This general method in implemented in
  the \pkg{NMF} package. In this case a single regression
  is fitted using the function \code{\link[NMF]{fcnnls}}.
}
\examples{
# random target matrix
x <- rmatrix(100, 20)
# random cell signatures
s <- rmatrix(100, 3)

# deconvolve using standard least-squares
res <- ged(x, s, 'lsfit')
coef(res)
# signatures are not updated
identical(basis(res), s)
\dontshow{
	stopifnot(identical(basis(res), s))
	stopifnot( nmf.equal(res, ged(x, s, 'lsfit')) )
}

# Fitting with fcnnls
res <- ged(x, s, 'lsfit', fit = 'nnls')
coef(res)
# signatures are not updated
identical(basis(res), s)
\dontshow{
	stopifnot(identical(basis(res), s))
	stopifnot( nmf.equal(res, ged(x, s, 'lsfit', fit = 'nnls')) )
}
}
\references{
  Abbas AR, Wolslegel K, Seshasayee D, Modrusan Z and Clark
  HF (2009). "Deconvolution of blood microarray data
  identifies cellular activation patterns in systemic lupus
  erythematosus." _PloS one_, *4*(7), pp. e6098. ISSN
  1932-6203, <URL:
  http://dx.doi.org/10.1371/journal.pone.0006098>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/19568420>.

  Van Benthem M and Keenan MR (2004). "Fast algorithm for
  the solution of large-scale non-negativity-constrained
  least squares problems." _Journal of Chemometrics_,
  *18*(10), pp. 441-450. ISSN 0886-9383, <URL:
  http://dx.doi.org/10.1002/cem.889>, <URL:
  http://doi.wiley.com/10.1002/cem.889>.

  Kim H and Park H (2007). "Sparse non-negative matrix
  factorizations via alternating non-negativity-constrained
  least squares for microarray data analysis."
  _Bioinformatics (Oxford, England)_, *23*(12), pp.
  1495-502. ISSN 1460-2059, <URL:
  http://dx.doi.org/10.1093/bioinformatics/btm134>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/17483501>.
}

