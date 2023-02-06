# DataSource definition for fetching data from CellMix registry
# 
# Author: Renaud Gaujoux
# Created: 08 May 2013
###############################################################################

#' @include DataSource.R
NULL

#' @S3method DataSource CellMix
DataSource.CellMix <- function(x, ...){
    list(get = function(key, destdir = NULL, ...){
                # NB: destdir is not used here
                ExpressionMix(key, ...)
            }
            , list = function(all = FALSE){ gedData(all = all) }
            , name = class(x)[1L]
    )
}
