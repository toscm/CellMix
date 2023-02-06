# Functions to handle Datasets
# 
# Author: Renaud Gaujoux
# Created: Apr 1, 2013
###############################################################################

ldata <- function(..., package='CellMix'){
    pkgmaker::ldata(..., package=package)
}

#' @include DataSource-ArrayExpress.R
#' @include GEOquery.R
NULL

#' Downloading Gene Expression Deconvolution Datasets
#'  
#' \code{GEDdownload} downloads gene expression datasets that 
#' for datasets that have an entry in the \pkg{CellMix} package's
#' internal registry. 
#' 
#' @param x access key of a registered dataset
#' @param ... extra arguments passed to \code{\link{download.file}} and/or
#' \code{\link{getGSE}}.
#' @param destdir destination directory for all the downloaded files
#' @param datasource Data source from where to fetch the dataset, e.g., 
#' an online dataset repository lke \emph{GEO} or \emph{ArrayExpress}.
#' If missing, it is determined automatically from the access key, but is honoured if provided.
#' Currently only values \code{'GEO'} and \code{'ArrayExpress'} (or its alias \code{'AE'}) are supported.
#' @param annotation annotation package (a string) to attach to the
#' created ExpressionSet object. 
#' @param verbose logical or numeric that indicate the verbosity level
#' 
#' @return the value returned by \code{\link{getGSE}}. 
#' 
#' @export
#' 
GEDdownload <- function(x, ..., destdir=GEDtmp(), datasource=NULL, annotation=NULL, verbose=TRUE){
    
    if( !is.verbose() || !missing(verbose) ){
        ol <- lverbose(verbose)
        on.exit( lverbose(ol) )
    }
    verbose <- lverbose()
    
    # get registry entry
    dobj <- gedData(x, error=FALSE)
    
    if( is.GEDdata(dobj) ){
        x <- dobj
        if( grepl("^(https?)|(ftp)", x$url) ){
            locfile <- file.path(destdir, basename(x$url))
            if( !file.exists(locfile) ){
                download.file(x$url, destfile = locfile, ...)
            }else vmessage("# Using cached version: '", locfile, "'")
            vmessage("# Loading file '", locfile, "'")
            return( getGSE(filename = locfile, destdir=destdir, ..., annotation = annotation) )
        }
        key <- x$key
    }else if( isString(x) && is.null(dobj) ){# not found in the local registry: try online databases
        key <- x
    }else{
        stop("Could not download dataset from the given specification [", class(x), ']')
    }
    
    # Fetch data from datasource 
    datasource <- DataSource(datasource, key = key)
    vmessage("# Fetching dataset '", key, "' from ", datasource$name)
    vmessage("# Local cache directory: '", destdir, "'")
    eset <- datasource$get(key, destdir = destdir, ...)
    
    # set annotation if provided
    if( !is.null(annotation) && is(eset, 'ExpressionSet') ){ 
        annotation(eset) <- annotation
    }
    
    # return ExpressionSet object
    eset
} 


# Returns the url associated with a given GED data.
GEDurl <- function(x){
    
    # get registry entry
    x <- gedData(x)
    
    if( grepl("^((https?)|(ftp))://", x$url) )
        return(x$url)
    key <- x$key
    # determine data source from accession key
    ds <- DataSource(key = key)
    
    if( ds$name == 'GEO' ){# GEO accession number
        paste0('http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', key) 
    }else if( ds$name == 'ArrayExpress' ){
        paste0("http://www.ebi.ac.uk/arrayexpress/experiments/", key)
    }else
        stop("Could not identify data provider for accession key '", key, "'.")
}

#' Internal Functions for CellMix Data Registry 
#' 
#' These functions are used internally to work with the gene expression dataset
#' registry within the CellMix package.
#' 
#' \code{GEDsave} saves the primary dataset loaded in the registry on disk. 
#' 
#' @param x GEDdata_entry object or list of GEDdata_entry objects.
#' If missing then all entries are used.  
#' @param force logical to force saving dataset(s)
#' 
#' @rdname GEDsave
#' @keywords internal
GEDsave <- function(x, force = FALSE){
    
    .save <- function(d_entry){
        #print(d_entry)
        p <- esetPath(d_entry)
        if( (force || !file.exists(p)) && !is.null(es <- eset(d_entry, load = FALSE)) ){
            message("Saving expression data from ", d_entry$key ," in cache file '", p, "'")
            saveRDS(es, p)
        }
    }
    
    # save all if no data is passed
    if( missing(x) ) x <- gedDataInfo(FALSE)$get_entries()
    
    if( is.GEDdata(x) ) .save(x)
    else if( is.list(x) ) lapply(x, GEDsave, force=force)
    else stop("invalid argument `x`: must be a list or a GEDdata_entry object.")
    
    invisible()
}

setOldClass('GEDdata_entry')

#' User Data Directory
#'  
#' \code{userData} returns the path to a local directory where package-related user data can be stored.
#' Note that a directory is \strong{always} created if necessary (see details).
#' 
#' If in interactive mode, the user is asked if the directory can be created in his home directory,
#' otherwise, or if the user does not allow the creation in his home, the directory is created 
#' in the current R session's temporary directory.  
#' 
#' @param ... path parts passed to \code{\link{file.path}} to be appended to 
#' the main path.
#' @param create logical that indicates if the directory should be created if it does not exists.
#' @param package name of the package associated with the user data path.
#' It is used to prefix the path, within the user R data directory. 
#' 
#' @seealso \code{\link{tempdir}}
#' @export
userData <- local({
            function(..., create=TRUE, package='base'){
                
                p <- file.path(Sys.getenv('HOME'), 'R-data', package)
                
                # ask the user about creating the directory
                if( create && !file.exists(p) ){
                    ans <- askUser(str_c("The ", package, " user data directory '", p, "' doen't exist. Do you want to create it?")
                            , idefault='y', default='y')
                    if( ans == 'n' ){
                        p <- file.path(tempdir(), 'R-data', package)
                    }
                    if( !file.exists(p) ){
                        message("Creating user data directory '", p, "'")
                        dir.create(p, recursive=TRUE)
                    }
                }
                file.path(p, ...)
            }
        })

#' \code{GEDpath} returns/sets the path to the local directory where CellMix data (e.g., cache) are stored.
#' 
#' @inheritParams userData
#' @param reset new value for the local directory
#' 
#' @rdname GEDsave 
#' @export
GEDpath <- local({
            .registryPath <- NULL
            function(..., create=TRUE, reset=NULL){
                # initialise path at first call
                if( is.null(.registryPath) || !is.null(reset) ){
                    p <- if( is.character(reset) ) reset[1L] 
                            else userData(create=create, package='CellMix')
                    .registryPath <<- p
                }
                file.path(.registryPath, ...)
            }
        })

#' \code{GEDtmp} returns the path to the local directory where downloaded data from GEO 
#' are stored.
#' 
#' @rdname GEDsave
GEDtmp <- function(..., create=TRUE){
    p <- GEDpath('GEO', ...)
    if( !file.exists(p) && create ){
        dir.create(p)
    }
    p
}

#' \code{GEDcache} adds an ExpressionSet object or an rds file to the local cache.
#' 
#' @param key unique dataset identifier. 
#' For datasets that are registered in the \pkg{CellMix} data registry, this  
#' correspond to the entry key.
#' @param object \code{ExpressionSet} object to add to cache or the path to
#' an existing RDS file that contains such an object. 
#' 
#' @export
#' @rdname GEDsave
GEDcache <- function(key, object){
    
    dest <- esetPath(key)
    # copy file
    if( isString(object) ){
        if( !is.file(object) ){
            stop("Could not copy file to GED cache: file '", object, "' does not exist")
        }
        file.copy(object, dest, overwrite = TRUE)
    }else{ # save object
        saveRDS(object, file=dest)
    }
}
