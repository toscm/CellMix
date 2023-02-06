# Data sources
# 
# Author: Renaud Gaujoux
# Created: 08 May 2013
###############################################################################

#' Gene Expression Data Sources
#' 
#' @description 
#' \code{DataSource} is an S3 generic function that creates prototype-like
#' objects whose \code{get} method fetches datasets from a given source/repository.
#' 
#' The \pkg{CellMix} package provides DataSource objects for GEO
#' (via the \pkg{GEOquery} package), ArrayExpress (via the \pkg{ArrayExpress} 
#' package), as well as the internal CellMix registry -- which itself fetches datasets 
#' from multiple repositories using other \code{DataSource} objects.
#' 
#' @param x data source key or object
#' @param ... extra argument to allow extension and passed to 
#' @export
DataSource <- function(x, ...){
    UseMethod('DataSource')
}

#' \code{isDataSource} tests if an object is a \code{DataSource}.
#' 
#' @export
#' @rdname DataSource
isDataSource <- function(x){
    is(x, 'DataSource')
}

asDataSource <- function(x, name, ...){
    
    # check type
    if( !is.list(x) ) stop("Invalid object type: DataSource objects are lists [", class(x), ']')
    
    .self <- x
    S3Class(.self) <- c(name, 'DataSource')
    
    # add has method if not present and a list method exists
    if( is.null(.self$has) && !is.null(.self$list) ){
        .self$has <- function(y, ...){
            y %in% .self$list(...)
        }
    }
    
    # return object
    .self
}

#' @S3method DataSource default
DataSource.default <- function(x = c('GEO', 'ArrayExpress', 'AE', 'CellMix'), key = NULL, quiet=FALSE, ...){
    
    # aliases
    aliases <- c(AE = 'ArrayExpress')
    if( was_missing <- missing(x) ) x <- NULL
    
    all_ds <- setdiff(gsub("DataSource\\.", "", methods('DataSource')), c('default', 'GEDdata_entry', 'DataSource'))
    
    # check validity and aliases for datasource
    if( !is.null(x) ){
        
        # return a list of DataSource objects
        if( length(x) > 1L){
            res <- sapply(x, DataSource, key = key, quiet = quiet, simplify = FALSE)
            S3Class(res) <- 'DataSourceList'
            return(res)
        }
        
        x <- match.arg(x, c(all_ds, names(aliases)))
        
        # resolve aliases
        if( x %in% names(aliases) ){
            x <- setNames(aliases[x], NULL)
        }
    }else{
        if( !is.null(key) ) x <- NULL # will be infered by key
        else{ # list all data sources
            if( !was_missing ) return()
            return(all_ds)
        }
    }
    
    # infer from key if datasource not already defined at this stage 
    if( is.null(x) && !is.null(key) ){
        x <- 
                if( grepl("^GSE", key) ) 'GEO' # GEO accession number			
                else if( grepl("^E-((MTAB)|(TABM)|(GEOD)|(MEXP))-", key) ) 'ArrayExpress' # ArrayExpress accession number
                else if( !quiet ){
                    stop("Could not identify data provider of accession key '", key, "'.")
                }else return()
    }
    
    # TODO: define S3 structure for data sources (name, regexp key pattern, url prefix, download method, etc...)
    # return data source
    S3Class(x) <- x
    asDataSource(DataSource(x), name = x)
}

#' @S3method DataSource DataSource
DataSource.DataSource <- function(x, ...) x

#' @S3method print DataSource
print.DataSource <- function(x, ...){
    cat("<DataSource: ", x$name,">\n", sep='')
    s <- capture.output(args(x$get))
    cat("  $get = ", paste0(strwrap(paste0(s[-length(s)], collapse = ''), exdent = 10), collapse = '\n'), "\n", sep='')
}

#' Fetching Data from Data Sources
#' 
#' \code{fetchData} retrieves some data from a given data source.
#' 
#' @param x dataset accession number (e.g. \code{'GSE12345'})
#' @param ... extra arguments passed to the datasource's \code{get} method.  
#' @param datasource datasource where the data is fetched from.
#' If \code{NULL}, then all known data sources are tried
#' 
#' @export 
fetchData <- function(x, ..., datasource = NULL){
    
    # get datasource definition(s)
    if( !isDataSource(datasource) ){
        datasource <- DataSource(datasource, key = x)
    }
    
    if( is(datasource, 'DataSourceList') ){ # try sequentially each data source
        for(d in datasource){
            res <- fetchData(x, d, ...)
            if( !is.null(res) ){
                return(res)
            }
        }
    }else{
        if( !is.null(datasource$has) ){
            if( datasource$has(x) ) datasource$get(x, ...)
        }else{
            res <- try(datasource$get(x, ...), silent = TRUE)
            if( is(res, 'try-error') ) return()
            return( res )
        }
    }
    
}

