# Searching identifiers
# 
# Author: Renaud Gaujoux
# Created: May 1, 2013
###############################################################################


#' @include Bioc.R
NULL

searchIDs <- function(ids, to, from, ...){
    map <- biocann_map(from, to)
    ids2 <- grep(ids, keys(map[[1]]), value=TRUE)
    convertIDs(ids2, to, from, ...)
}
