# DataSource definition for fetching data from ArrayExpress
# 
# Author: Renaud Gaujoux
# Created: 08 May 2013
###############################################################################

#' @include DataSource.R
NULL

#' @S3method DataSource ArrayExpress 
DataSource.ArrayExpress <- function(x, ...){
    list(get = getArrayExpress, name = class(x)[1L])
}

getArrayExpress <- function(key, destdir = GEDtmp(), cache = TRUE, ..., verbose = FALSE){
    
    # set context verbosity level if necessary
    if( missing(verbose) && is.verbose() ) verbose <- lverbose()
    
    # require ArrayExpress
    uq_requirePackage('ArrayExpress', load=TRUE, msg=str_c("for downloading datasets from ArrayExpress."), ptype='BioCsoft')
    
    # try local
    if( lverbose() <= 1 ){ # show warnings if high verbose level
        opt <- options(warn=-1L)
        on.exit( options(opt) )
    }
    vmessage("# Checking local files ...")
    
    # suppress output from getAE
    getAE <- silenceF(ArrayExpress::getAE, verbose)
    
    ae <- NULL
    if( cache ){# lookup local cache
        ae <- try(getAE(key, path = destdir, type='processed', local=TRUE), silent=TRUE)
        if( is(ae, 'try-error') ){
            ae <- NULL
            vmessage("# FAILED")
        }else vmessage("# OK")
    }
    
    if( is.null(ae) ){ # get from remote repo
        vmessage("# Downloading data files ...")
        ae <- getAE(key, path=destdir, type='processed', local=FALSE)
        vmessage("# DONE")
    }
    # FIXME: ArrayExpress - procset crashes if multiple processed files are present, 
    # even if only one of them is actually needed.
    # Here we remove unnecessary files from the list
    # read sample annotation file
    if( (nfiles <- length(ae$processedFiles)) > 1 ){
        vmessage("# Limiting to necessary expression data files ... ", appendLF=FALSE)
        if( file.exists(sannot <- file.path(destdir, paste0(key, '.sdrf.txt'))) ){
            la <- readLines(sannot)
            fused <- sapply(ae$processedFiles, function(f) sum(grepl(f, la[-1], fixed=TRUE)) > 1)
            if( length(not_used <- which(!fused)) ){
                ae$processedFiles <- ae$processedFiles[-not_used]
            }
            vmessage("OK [", length(ae$processedFiles), "/", nfiles,"]")
        }else vmessage("NO [no sample annotation]")
    }
    #
    vmessage("# Loading/Processing dataset ...")
    cnames = getcolproc(ae)
    #		print(ae)
    #		print(cnames)
    procset <- silenceF(ArrayExpress::procset, verbose)
    eset <- procset(ae, cnames[2])
    vmessage("# DONE")
    
    # FIXME: ArrayExpress - procset somteimes adds an extra column full of NAs
    if( all(is.na(exprs(eset)[, ncol(eset)])) ){
        vmessage("# Dropping last invalid column ...")
        eset <- eset[, -ncol(eset)]
        vmessage("OK [", ncol(eset), ']')
    }
    # FIXME: ArrayExpress makes a bad job in loading the sample annotation data
    if( !ncol(pData(eset)) ){
        vmessage("# Fixing sample phenotypic annotation ...")
        readPhenoData <- silenceF(ArrayExpress:::readPhenoData, verbose)
        pd <- readPhenoData(ae$sdrf, destdir)
        pd <- pData(pd)
        sid <- sampleNames(eset)
        if( is_NA(idc <- matchColumn(sid, pd)) ){
            if( is_Illumina(featureNames(eset)) ){ # try removing X prefix
                sid <- substring(sid, 2)
                idc <- matchColumn(sid, pd)
            }
        }
        if( !is_NA(idc) ){
            i <- match(sid, pd[[idc]])
            if( !any(is.na(i)) ){
                pd <- pd[i, ]
                rownames(pd) <- sid
                pData(eset) <- pd
                vmessage("OK [attached ", length(varLabels(eset)), ' variables]')	
            }else{
                vmessage("FAILED [missing some sample annotations]\n  NOTE: Check that the annotation file '", ae$sdrf, "' contains data for all samples.")
            }
        }else vmessage("FAILED [could not match samples]\n  NOTE: Check that the annotation file '", ae$sdrf, "' contains data for all samples.")
    }
    #
}
