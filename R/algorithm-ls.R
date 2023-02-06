# Deconvolution algorithms based on standard least-squares:
#
# Simple strategy from Abbas et al. (2009)
# 
# Author: Renaud Gaujoux
# Creation: 18 Jan 2012
###############################################################################

#' @include registry-algorithms.R
NULL

#' Partial Gene Expression Deconvolution by Least-Square 
#' 
#' 
#' \code{.gedLSfit} provides access to partial deconvolution methods that are 
#' based on least-squares fits. 
#' 
#' @inheritParams .qprog
#' @param rescale logical used when estimating proportions from signatures, 
#' that indicates if the esti,ated coefficients should be scaled to sum up to 
#' one (\code{TRUE}) or left as estimated by the linear regression (\code{FALSE}).
#' This scaling is performed after the coefficients have been forced to be 
#' nonnegative.
#' @param fit least-square fitting method: \code{ls} uses \code{\link{lm}}, 
#' \code{nnls} uses \code{\link[NMF]{fcnnls}}. 
#' @param ... extra arguments passed to fitting the methods \code{\link{.nn_lsfit}} 
#' or \code{\link[NMF]{.fcnnls}}.
#' 
#' @return an \code{\link{NMF}} object.
#' 
#' @rdname gedLSfit
#' @keywords internal
.gedLSfit <- function(X, seed, rescale=TRUE, fit = c('ls', 'nnls'), ...){
    
    fit <- match.arg(fit)
    if( hasBasis(seed) && !hasCoef(seed) ){
        
        vmessage("Estimating cell proportions from cell-specific signatures [lsfit: ", fit, "]")
        # fit with lsfit
        if( fit == 'ls' ){
            res <- .nn_lsfit(x = basis(seed), y = X, ...)
            # estimates
            res_coef <- res$coef
            # fit object
            res_fit <- res$fit
        }else if( fit == 'nnls' ){
            res <- fcnnls(basis(seed), X, ...)
            # estimates
            res_coef <- res$x
            # fit object
            res_fit <- res
        }
        # store fit
        seed$coeffit <- res_fit
        # update coef
        coef(seed) <- res_coef
        # rescale coef
        if( rescale ) coef(seed) <- scoef(seed)
        
    }else if( !hasBasis(seed) && hasCoef(seed) ){
        
        vmessage("Estimating cell-specific signatures from cell proportions [cs-lsfit: ", fit, "]")
        
        if( fit == 'ls' ){
            # fit the transposed problem with lsfit
            res <- .nn_lsfit(x=t(coef(seed)), y=t(X), ...)
            # estimates
            res_coef <- t(res$coef)
            # fit object
            res_fit <- res$fit
        }else if( fit == 'nnls' ){
            # fit the transposed problem with fcnnls
            res <- fcnnls(t(coef(seed)), t(X), ...)
            # estimates
            res_coef <- t(res$x)
            # fit object
            res_fit <- res
        }
        # store fit
        seed$basisfit <- res_fit
        # update basis 
        basis(seed) <- res_coef
        
    }else{
        stop("lsfit - Only partial fitting is currently supported (i.e. either the signatures or the proportions must be provided).")
    }
    
    # return updated seed
    seed
}

#' \code{.nn_lsfit} implements a standard least-square fit with various 
#' different procedures to enforce nonnegative coefficients.
#' In particular, it implements the iterative procedure described in 
#' \cite{Abbas2009}. 
#' 
#' @param x matrix of known cell-specific profiles (i.e. cell signatures), 
#' with features in rows and cell type in columns.
#' @param y matrix of observed mixed expression data, with features in rows and
#' samples in columns.
#' The number of samples must be greater -- or equal -- than the number of 
#' cell types.  
#' @param nneg specification of the method used to enforce the nonnegativity 
#' of the estimated proportions.
#' 
#' Accepted values are:
#' \describe{
#' \item{\code{'iterate'}:}{ applies the procedure described in \cite{Abbas2009}. 
#' For each sample separately, a sequence of least-square fits are performed, 
#' starting with all cell types, and where the cell type corresponding to the 
#' lowest negative fitted coefficient is excluded from the next fit, and its 
#' associated final proportion set to zero.
#' This iterative process stops when all coefficients are nonnegative.}
#' \item{\code{'pmax'}:}{ single least-square fit, where all negative 
#' estimated proportions are set to zero.}
#' \item{\code{NA} or \code{'none'}: }{ single least-square fit where 
#' the estimated proportions are returned unconstrained.}
#' } 
#' 
#' @rdname gedLSfit
.nn_lsfit <- function(x, y, nneg=c('iterate', 'pmax', 'none'), ...){    
    
    # by default enforce nonnegativity only if the input data contains negative values
    if( missing(nneg) && min(y) < 0 ) nneg <- NA
    if( !is_NA(nneg) ){
        nneg <- match.arg(nneg)
    }
    
    # fit with lsfit
    fit0 <- lm(y ~ -1 + x, ...)
    # store initial fit
    res <- list(fit = fit0)
    
    # enforce nonnegative values if necessary
    if( !is_NA(nneg) ){
        if( nneg == 'pmax' ) res$coef <- pmax(coef(fit0)) 
        else if( nneg == 'iterate' ){
            r <- ncol(x)
            # iterate as described in Abbas et al. (2009)
            res$coef <- 
                sapply(seq(ncol(coef(fit0))), function(k, ...){
                            a <- coef(fit0)[, k]
                            i <- integer()
                            while( any(a<0) && length(i) < r ){
                                i <- c(i, which.min(a))
                                a[i] <- 0
                                fit <- lm(y[,k] ~ -1 + x[, -i, drop=FALSE], ...)
                                stopifnot( length(a[-i]) == length(coef(fit)) )
                                a[-i] <- coef(fit)
                            }
                            a
                        })
        }
    }
    
    # re-apply dimnames
    rownames(res$coef) <- colnames(x)
    colnames(res$coef) <- colnames(y)
    
    res
}



# Registration of NMF method 'lsfit'
nmfAlgorithm.lsfit <- setNMFMethod('lsfit', .gedLSfit
        , objective='euclidean' 
        , mixed=TRUE
        , overwrite=TRUE
)

#' Partial Gene Expression Deconvolution by Standard Least-Squares
#' 
#' Estimates cell/tissue proportions given a known set of cell/tissue-specific 
#' expression signatures, using standard least-squares as proposed by 
#' \cite{Abbas2009}.
#'
#' The default algorithm uses a heuristic to enforce the nonnegativity of the estimated 
#' proportions, that consists in fitting successive regressions, each time excluding 
#' the most negative coefficient from the model, until all coefficients are nonnegative.
#' In this case all regressions are fitted using the function \code{\link{lm}}.
#' 
#' An alternative least-square fitting method is included for test/experimental 
#' purposes.
#' It uses the fast combinatorial nonnegative least-square method of 
#' \cite{VanBenthem2004}, which was adapted by \cite{KimH2007} to perform 
#' nonnegative matrix factorization of gene expression 
#' -- but not originally for deconvolution. 
#' This general method in implemented in the \pkg{NMF} package.
#' In this case a single regression is fitted using the function 
#' \code{\link[NMF]{fcnnls}}.
#'  
#' @inheritParams .gedLSfit
#' @aliases lsfit-ged
#' @cite Abbas2009
#' @examples 
#' 
#' # random target matrix
#' x <- rmatrix(100, 20)
#' # random cell signatures
#' s <- rmatrix(100, 3)
#' 
#' # deconvolve using standard least-squares
#' res <- ged(x, s, 'lsfit')
#' coef(res)
#' # signatures are not updated
#' identical(basis(res), s)
#' \dontshow{ 
#' 	stopifnot(identical(basis(res), s))
#'	stopifnot( nmf.equal(res, ged(x, s, 'lsfit')) ) 
#' }
#' 
#' # Fitting with fcnnls
#' res <- ged(x, s, 'lsfit', fit = 'nnls')
#' coef(res)
#' # signatures are not updated
#' identical(basis(res), s)
#' \dontshow{ 
#' 	stopifnot(identical(basis(res), s))
#'	stopifnot( nmf.equal(res, ged(x, s, 'lsfit', fit = 'nnls')) ) 
#' }
#' 
gedAlgorithm.lsfit <- setGEDMethod(key='lsfit'
        , description = "Partial deconvolution of proportions using least-squares fits"
        , algorithm = 'lsfit' 
        , reqBasis = TRUE, outBasis = FALSE
        , reqCoef= FALSE, outCoef = TRUE
        , reqMarker = FALSE
        , maxIter = 1L
        , cite = "Abbas2009"
)

#' Cell-Specific Expression by Standard Least-Squares
#' 
#' Estimates cell-specific proportions given known proportions 
#' expression signatures, using least-squares fitting.
#'
#' The algorithm applies the same methods as the ged algorithm 
#' \code{\link[=lsfit-ged]{lsfit}} but to the transposed problem of 
#' estimating signatures from proportions.
#' It is included in the \pkg{CellMix} package for test/experimental purposes. 
#'  
#' @inheritParams gedAlgorithm.lsfit
#' @aliases cs-lsfit-ged
#' @examples 
#' 
#' # random target matrix
#' x <- rmatrix(100, 20)
#' # random cell proprtions
#' p <- rmatrix(3, 20)
#' 
#' # deconvolve using standard least-squares
#' res <- ged(x, p, 'cs-lsfit')
#' head(basis(res))
#' # used proportions are stored in the result coefficient matrix
#' identical(coef(res), p)
#' 
#' # deconvolve using nonnegative least-squares
#' res <- ged(x, p, 'cs-lsfit', fit = 'nnls')
#' head(basis(res))
#' # used proportions are stored in the result coefficient matrix
#' identical(coef(res), p)
#' 
gedAlgorithm.cs_lsfit <- setGEDMethod(key='cs-lsfit'
        , description = "Partial deconvolution of cell signatures using least-squares fits"
        , algorithm = 'lsfit'
        , reqBasis = FALSE, outBasis = TRUE
        , reqCoef= TRUE, outCoef = FALSE
        , reqMarker = FALSE
        , maxIter = 1L
)
