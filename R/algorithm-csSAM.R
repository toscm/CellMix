# csSAM wrapper
# 
# Wraps call to the csSAM package from Shen-Orr et al. (2010)
# Author: Renaud Gaujoux
###############################################################################

#' @include registry-algorithms.R
#' @include plots.R
NULL

# simple warning note
wnote <- function(..., immediate. = TRUE){
    warning(..., immediate.=immediate., call.=FALSE)
}

#' Cell-specific Differential Expression with csSAM
#' 
#' This function is adapted from the function \code{\link[csSAM]{csSamWrapper}} 
#' in the \code{\link{csSAM}} package, to integrate the csSAM algorithm from 
#' \cite{Shen-Orr2010} into the \pkg{CellMix} framework of deconvolution algorithms.
#' 
#' @inheritParams csSAM::csSamWrapper
#' @param Y target global gene expression matrix (n x p), with samples in columns, ordered in the 
#' same order at the cell proportions data in \var{x}.
#' @param x known cell proportions as a matrix (k x p) or an \code{\linkS4class{NMF}} model 
#' containing the cell proportions in the coefficient matrix -- and a normally
#' empty basis matrix.
#' The proportions must be ordered in the same order as the samples in the target matrix.
#' 
#' For \code{csTopTable}, a csSAM fit as return by \code{\link{ged}}.
#' @param data specification of the sample groups.
#' If not missing, it must be a factor or coercible to a factor, with length the 
#' number of samples, i.e. columns, in the target matrix. 
#' @param nperms The number of permutations to perform.
#' It is only used when computing cell-specific differential expression between
#' groups specified in argument \code{data}.  
#' @param verbose logical that indicates if verbose messages should be shown.
#' 
#' @return Returns an NMF model object, with the cell specific differential 
#' expression stored in the \code{\link{basis}} matrix.
#' 
#' The following details about the fit can be extracted using \code{\link{basisfit}}:
#' \item{csfit}{A list object containing a fit (cell-type specific differential
#' expression) for each sample group. 
#' Each element in the list is an object returned by \code{\link[csSAM]{csfit}}.}
#' \item{csSAMfit}{ A list output of the fdrCsSAM function, containing data 
#' from the FDR computation. These are used by the \code{\link{csplot}} and 
#' \code{\link{csTopTable}}} methods associated with \code{csSAM} fits.
#' 
#' @author 
#' Original function: Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' 
#' Adaptation for \pkg{CellMix}: Renaud Gaujoux
#' @seealso
#' \code{\link{csfit}},\code{\link{csSAM}},\code{\link{fdrCsSAM}}
#' @cite Shen-Orr2010
#' @import csSAM
.csSAM <- function(Y, x, data=NULL
        , nperms = 200, alternative = c('all', 'two.sided', 'greater', 'less')
        , standardize=TRUE, medianCenter=TRUE
        , logRm =FALSE, logBase = 2, nonNeg=TRUE
        , verbose = lverbose()){ 
    #function(G,cc,y,nperms = 200,alternative = 'two.sided'
#				,standardize=TRUE,medianCenter=TRUE, logRm =FALSE,logBase = 2,nonNeg=TRUE
#				,fileName='csSAMout.pdf') {
    
    ## map arguments to those of csSamWrapper
    if( is.matrix(x) ){ # convert into an NMF model if necessary
        x <- nmfModel(Y, H=x, force.dim=FALSE)
    }
    
    # check and remove absent cell types
    if( any(nc <- apply(coef(x), 1L, function(x) all(x==0))) ){
        cn <- basisnames(x)[nc]
        wnote("Removed absent cell types from the estimation: "
                , str_out(cn), ' [', length(cn), ']')
        # subset on basis
        x <- x[!nc]
    }
    #
    
    if( isExpressionSet(Y) ) Y <- exprs(Y)
    # data contains the sample group definition
    if( is.null(data) ){
        data <- factor(rep(1, ncol(Y)))
    }else{
        if( !is.factor(data) ) data <- factor(data, levels=unique(data))
        # remove absent levels
        data <- droplevels(data)
    }
    # transpose/copy for csSAM
    G <- t(Y); cc <- t(coef(x)); y <- as.integer(data)
    ##
    
    # remove NA-labeled samples
    if( length(i_rm <- which(is.na(y))) ){
        wnote('Dropped NA-labeled samples: ', str_out(rownames(G)[i_rm]), ' [', length(i_rm), ']')
        y <- y[-i_rm]
        cc <- cc[-i_rm, , drop=FALSE]
        G <- G[-i_rm, , drop=FALSE]
    }
    
    deconv <- list()
    
#	numset = length(unique(y))
    n <- summary(data, maxsum=Inf)
    numset <- nlevels(data)
    
    vmessage("Groups: ", str_out(n, Inf, use.names=TRUE, sep=' | '))
    # check parameters
    if( numset > 2L ){
        stop("csSAM - Cannot handle more than 2 groups of samples [", numset, "]")
    }
    if( nrow(G) != nrow(cc) ){
        stop("csSAM - Incompatible dimensions: number of cell proportions [", nrow(cc),"]"
                , " should match the number of samples [", nrow(G),"]")
    }
    #
    
    numgene = ncol(G)
    numcell = ncol(cc)
    cellNames = colnames(cc)
    
#	n <- vector(mode = "logical", length  = numset)
#	for (i in 1:numset) n[i] =sum(y==i)
    
#	for (curset in 1:numset) {
#		deconv[[curset]]= csfit(cc[y==curset,], G[y==curset,],logRm,logBase)
#	}
    vmessage("Fitting cell-specific linear model ... ", appendLF=FALSE)
    sets <- split(1:nrow(G), y)
    deconv <- lapply(sets, function(i){
                csfit(cc[i,,drop=FALSE], G[i,,drop=FALSE], logRm = logRm, logBase = logBase)
            })
    vmessage('OK')
    
    # early exit if only one group
    if( numset == 1L ){
        d <- dimnames(x)
        basis(x) <- t(deconv[[1L]]$ghat) 
        dimnames(x) <- d
        x$residuals <- t(deconv[[1L]]$residuals)
        x$se <- t(deconv[[1L]]$se)
        x$rawfit <- deconv[[1L]]
        return(x)
    }
    
    #rhat <- array(dim = c(numcell,numgene))
    vmessage("Computing csSAM model statistics ... ", appendLF=FALSE)
    rhat <- csSAM(deconv[[1]]$ghat, deconv[[1]]$se,
            n[1], deconv[[2]]$ghat, deconv[[2]]$se, n[2],
            standardize, medianCenter, nonNeg)
    vmessage("OK")
#	tt.sam <- runSAM(G, y)
    alternative <- match.arg(alternative)
    vmessage("Computing fdr using ", nperms, " permutations ... ", appendLF=FALSE)
    elapsed <- system.time({
                fdr.csSAM <- fdrCsSAM(G,cc,y,n,numcell,numgene, rhat
                        , nperms = nperms, alternative = alternative
                        , standardize = standardize, medianCenter = medianCenter
                        , logRm = logRm, logBase = logBase, nonNeg = nonNeg
                        , verbose = max(verbose - 1L, 0L))
                fdr.csSAM$alternative <- alternative
            })
    vmessage("OK")
    lmessage(2L, "Timing:")
    if( verbose >= 2L ) show(elapsed)
    
#	fdr.sam <- fdrSAM(G, y, nperms=nperms, tt.sam, alternative)
#    vmessage("Finding signature genes ... ", appendLF=FALSE)
#    sigGene <- findSigGene(G, cc, y, rhat, fdr.csSAM)
#    vmessage('OK')
    
#	plotCsSAM(fdr.csSAM, fdr.sam,alternative,cellNames,numcell, file = fileName)
#	return(list(deconv = deconv, fdr.csSAM = fdr.csSAM, fdr.SAM = fdr.sam,sigGene.csSAM = sigGene,fileName = fileName))
    
    # wrap the result into the NMF model
    # remove permuation results if wrapping (to avoid memory issue)
    if( .cellmix.getOption('wrap') ) fdr.csSAM$rhatperm <- NULL
    x$basisfit <- list(csfit = deconv, csSAMfit = fdr.csSAM
#                    , sigGene.csSAM = sigGene
#					, fdr.SAM = fdr.sam
    )
    b <- t(rhat)
    rownames(b) <- rownames(Y)
    colnames(b) <- basisnames(x)
    basis(x) <- b
    x$call <- match.call()
    x$nperms <- nperms
    x$data <- data
    
    # return updated NMF model
    x
}

fdrCsSAM <- function (G,cc,y,n,numcell,numgene,rhat,nperms,alternative='two.sided',standardize=TRUE,medianCenter=TRUE,logRm=FALSE,logBase = 2,nonNeg=FALSE, verbose = TRUE) {
    numgene=ncol(G)
    rhatperm <- array(dim = c(nperms,numcell,numgene))
    progress <- iterCount(n=nperms, title='', verbose = verbose)
    lapply(1:nperms, function(i){
        progress()
        ystar = y[sample.int(length(y))]
        perm = list()
        for (curset in 1:2) {
            perm[[curset]]= csfit(cc[ystar==curset,], G[ystar==curset,],logRm,logBase)
        }
        rhatperm[i,,] <<- csSAM(perm[[1]]$ghat, perm[[1]]$se, n[1], perm[[2]]$ghat, perm[[2]]$se, n[2],
                standardize, medianCenter, nonNeg)
        NULL
    })
    progress(nperms, appendLF = FALSE)
    
    # compute fdr for all alternatives
    FDR <- list()
    if( alternative != 'all' ){
        FDR[[alternative]] <- csSAM_fdr(rhat, rhatperm, alternative)
    }else{
        FDR[['two.sided']] <- csSAM_fdr(rhat, rhatperm, 'two.sided')
        FDR[['greater']] <- csSAM_fdr(rhat, rhatperm, 'greater')
        FDR[['less']] <- csSAM_fdr(rhat, rhatperm, 'less')
    }
    
    return (list(FDR = FDR, rhat = rhat, rhatperm = rhatperm, alternative = alternative))
}

csSAM_fdr <- function(rhat, rhatperm, alternative, numcell = nrow(rhat), nperms = dim(rhatperm)[1L]){
    
    vmessage("Alternative '", alternative,"' ... ", appendLF=FALSE)
    cutp.g=matrix(NA,nrow=numcell,ncol=100)
    numcut = ncol(cutp.g)
    
    fdr.g=ncall.g=nperm.g<-array(dim = c(numcell, numcut))
    
    for(j in 1:numcell)
        cutp.g[j,]=seq(0,max(abs(rhat[j,])),length=100)	
    for (i in 1:numcut) {
        for (curcell in 1:numcell) {
            if(alternative == 'two.sided') {
                fdr.g[curcell,i]=sum(abs(rhatperm[,curcell,])>cutp.g[curcell,i])/nperms /sum(abs(rhat[curcell,])>cutp.g[curcell,i])
                ncall.g[curcell,i]=sum(abs(rhat[curcell,])>cutp.g[curcell,i])
            }
            if(alternative == 'greater') {
                fdr.g[curcell,i]=sum(rhatperm[,curcell,]>cutp.g[curcell,i])/nperms /sum(rhat[curcell,]>cutp.g[curcell,i])
                ncall.g[curcell,i]=sum(rhat[curcell,]>cutp.g[curcell,i])
            }
            if(alternative == 'less') {
                # [RG] BUG FIX: should be < - cutp.g[curcell,i]
#			fdr.g[curcell,i]=sum(rhatperm[,curcell,]< -cutp.g[curcell,i])/nperms /sum(rhat[curcell,]>cutp.g[curcell,i])
                fdr.g[curcell,i]=sum(rhatperm[,curcell,]< -cutp.g[curcell,i])/nperms /sum(rhat[curcell,] < -cutp.g[curcell,i])
                ncall.g[curcell,i]=sum(rhat[curcell,]< -cutp.g[curcell,i])
            }		
        }
    }
    
    fdr.g =pmin(fdr.g,1)
    
    for (j in 1:numcell)	 {
        fdr.g[j,]=make.monotone(fdr.g[j,])
    }	
    
    vmessage('OK')
    list(fdr.g=fdr.g, cutp.g = cutp.g, ncall.g = ncall.g, alternative = alternative)
}

# compile functions
csSAM_fdr <- compiler::cmpfun(csSAM_fdr)
fdrCsSAM <- compiler::cmpfun(fdrCsSAM)

# overload csSAM function: use fast C++ code
# Will eventually be included in csSAM
#.findSigGene <- function (G, cc, y, rhat, csSAMData) 
#{
#	.Call('findSigGenes', rhat, csSAMData$cutp.g, csSAMData$fdr.g, PACKAGE='CellMix')
#}

# plain R version of findSigGene to test improvement in speed
.findSigGene_R <-
        function(G,cc,y,rhat,csSAMData) {
    numgene=ncol(G)
    numcell = ncol(cc)
    thresholdVec = csSAMData$fdr.g
    thresholdLen = length(thresholdVec[numcell,])
    sigGene <- array(dim = c(numcell, numgene))
    sigGene[,] = 1
    
    for (curThresh in 1:thresholdLen) {
        for (curcell in 1:numcell) {
            for (curgene in 1:numgene) {
                if(abs(rhat[curcell,curgene]) >= abs(csSAMData$cutp.g[curcell,curThresh])) {
                    sigGene[curcell,curgene] = csSAMData$fdr.g[curcell,curThresh]
                }
            }
        }
    }
    
    return (sigGene)
}

#' Partial Gene Expression Deconvolution with csSAM
#' 
#' Estimates cell/tissue proportions given a known set of cell/tissue-specific 
#' expression signatures, using standard least-squares, as implemented 
#' by the package \code{\link{csSAM}}.
#'
#' All regressions are fitted using the function \code{\link{lsfit}}.
#'  
#' @inheritParams .csSAM
#' @aliases csSAM-ged
#' @cite Shen-Orr2010
#' @examples 
#' 
#' # random global expression
#' x <- rmix(3, 100, 20)
#' basisnames(x) <- paste('Cell', 1:nbasis(x))
#' # extract true proportions
#' p <- coef(x)
#' 
#' # deconvolve using csSAM
#' res <- ged(x, p, 'csSAM')
#' head(basis(res))
#' # proportions are not updated
#' identical(coef(res), p)
#' \dontshow{ 
#' 	stopifnot(identical(coef(res), p))
#'	stopifnot( nmf.equal(res, ged(x, p, 'csSAM')) ) 
#' }
#' 
#' # estimate cell-specific differential expression between 2 groups
#' gr <- gl(2, 10)
#' res <- ged(x, p, 'csSAM', data = gr, nperms=20, verbose=TRUE)
#' head(basis(res))
#' # plot FDRs
#' csplot(res)
#' # extract fdr for top differentially expressed gene in each cell type
#' t <- csTopTable(res)
#' str(t)
#' 
gedAlgorithm.csSAM <- setGEDMethod(key='csSAM'
        , description = "Estimates cell/tissue specific signatures from known proportions using SAM"
        , algorithm = .csSAM
        , reqBasis = FALSE, outBasis = TRUE
        , reqCoef= TRUE, outCoef = FALSE
        , reqMarker = FALSE, outMarker = TRUE
        , maxIter = 1L
        , cite = "Shen-Orr2010" 
)

#' The S3 method \code{csTopTable} for csSAM fits returns, for each feature, the false discovery 
#' rates of differential expression between groups of samples within each cell type, 
#' as computed by \code{\link[=csSAM]{fdrCsSAM}} when running csSAM.
#' These are returned as a list, whith one element per cell type.
#' 
#' @seealso \code{\link[csSAM]{fdrCsSAM}}, \code{\link{csTopTable}}
#' 
#' @inheritParams csTopTable.matrix
#' 
#' @S3method csTopTable csSAM
#' @rdname gedAlgorithm.csSAM
#' 
csTopTable.csSAM <- function(x, alternative = c('two.sided', 'greater', 'less'), ...){
    
    # extract fitted model
    csSAMfit <- if( isNMFfit(x) ) basisfit(x)$csSAMfit else x
    
    if( is.null(csSAMfit$rhat) ){
        stop("Cannot compute top table: csSAM fit does not contain gene significance data.\n  Was csSAM run with a group variable?")
    }
    # vmessage("Finding signature genes ... ", appendLF=FALSE)
    alternative <- match.arg(alternative, names(csSAMfit$FDR))
    csSAMData <- csSAMfit$FDR[[alternative]]
    rhat <- csSAMfit$rhat
    G <- matrix(NA, 0, ncol(rhat))
    cc <- matrix(NA, 0, nrow(rhat))
    res <- findSigGene(G, cc, y = NULL, rhat, csSAMData)
    res <- t(res)
    # vmessage('OK')
    
    rownames(res) <- rownames(x)
    colnames(res) <- basisnames(x)
    
    csTopTable(res, ...)
}

#' The S3 method \code{csplot} for csSAM fits plots cell-specific fdr cumulative distributions. 
#' 
#' @param types index or names of the type to plot.
#' They need to be found in the fit data. 
#' @inheritParams graphics::plot.default 
#' 
#' @S3method csplot csSAM
#' @rdname gedAlgorithm.csSAM
#' @import ggplot2 
#' @importFrom plyr ldply
csplot.csSAM <- function(x, types=NULL, alternative = 'all', xlab='# called', ylab='FDR', ylim=c(0,1), ...){
    
    # extract fitted model
    csSAMfit <- if( isNMFfit(x) ) basisfit(x)$csSAMfit else x
    
    if( !is.list(csSAMfit) ) stop("Invalid input data: must be a list.")
    if( is.null(csSAMfit$rhat) ){
        warning("Nothing to plot: data does not contain any fdr data (did csSAM run with a group variable in argument `data`?)")
        return(invisible())
    }
    
    ncell <- nrow(csSAMfit$rhat)
    cellnames <- rownames(csSAMfit$rhat)
    
    if( is.null(cellnames) ){
        if( is.character(types) ){
            if( length(types) > ncell ){
                stop("Invalid number of types (", length(types), "):"
                        , " should have at most the same number of cell type data in `x` (", ncell, ")")
            }
            cellnames <- types
        }else cellnames <- as.character(seq(ncell))
    }
    if( is.null(types) ) types <- seq(ncell)
    if( is.character(types) ){ # partial match the types
        types <- pmatch(types, cellnames)
        types <- types[!is.na(types)]
        if( !length(types) )
            stop("None of the types ", str_out(types), " matched [available: ", str_out(cellnames),"]")
    }
    
#	pdf(file = fileName, height = 11, width = 8)
    op <- par(mfrow = mfrow(ncell))
    on.exit( par(op) )
    
#	plot(SAMdata$ncall.sam, SAMdata$fdr.sam, xlab = "## called", 
#			ylab = "FDR", type = "l", log = "x", ylim = c(0, 1))
#	title(paste("SAM", alternative))
    
    # extract alternative data
    if( alternative != 'all' ){
        alternative <- match.arg(alternative, names(csSAMfit$FDR))
        csSAMfit$FDR <- csSAMfit$FDR[alternative]
    }
    
    df <- ldply(types, function(i){                
        df <- ldply(csSAMfit$FDR, function(x){
                    j <- which(x$ncall.g[i, ]<=0)
                    data.frame(Alternative = x$alternative, x = rev(x$ncall.g[i, -j]), y = rev(x$fdr.g[i, -j]))
                })
        df[['Cell.type']] <- cellnames[i]
        df
    })

    # compute breaks in log scale
    breaks <- .log10_break(df$x)
    
    # create ggplot plot
    p <- ggplot(df, aes(x = x, y = y)) + 
            geom_line(aes(colour = Alternative)) +
            scale_y_continuous(ylab, limits=ylim) + 
            scale_x_continuous(xlab, breaks = breaks) +
            coord_trans(xtrans = "log2") +
            theme_bw() +
            facet_grid(~ Cell.type)
    return(p)
                
#    lapply(types, function(i){
#                j <- which(csSAMdata$ncall.g[i, ]<=0)
#                plot(csSAMdata$ncall.g[i, -j], csSAMdata$fdr.g[i, -j]
#                        , xlab = xlab, ylab = ylab
#                        , type = "l"
#                        , log = "x", ylim = ylim)
#                title(str_c(cellnames[i], ' - ', csSAMdata$alternative))
#            })
#    invisible(csSAMdata)
#	dev.off()
    
}

.log10_break <- function(x){
    m <- max(x)
    n <- 0
    while( m > 10 ){
        m <- ceiling(m/10)
        n <- n+1
    }
    r <- as.numeric(outer(c(1, m), 10^(0:n)))
    sort(r)
}
