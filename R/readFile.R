# Generic interface to load data from files/urls
# 
# Author: Renaud Gaujoux
# Created: 09 May 2013
###############################################################################


readFile <- function(file, ...){
    UseMethod('readFile')
}
readFile.default <- function(file, ...){
    S3Class(file) <- tolower(file_extension(file))
    readFile(file, ...)
}

writeFile <- function(file, ...){
    UseMethod('writeFile')
}
writeFile.default <- function(file, ...){
    S3Class(file) <- toupper(file_extension(file))
    writeFile(file, ...)
}

# RDS files
readFile.RDS <- function(file, ...) readRDS(file, ...)
writeFile.RDS <- function(file, ...) saveRDS(..., file = file)

# RDA
readFile.RDA <- function(file, ...){
    e <- new.env()
    load(file = file, envir = e)
    e
}
writeFile.RDA <- function(file, ...) save(..., file = file)

# CSV files
readFile.CSV <- function(file, ...) read.csv(file, ...)
writeFile.CSV <- function(file, ...) write.table(..., file = file, sep = ",")

# TSV files
readFile.TSV <- function(file, ...) read.delim(file, ..., sep = "\t")
writeFile.TSV <- function(file, ...) write.table(..., file = file, sep = "\t")
