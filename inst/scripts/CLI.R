# Command Line Interface for CellMix
# 
# Author: Renaud Gaujoux
# Created: May 7, 2013
###############################################################################

source('CLI-utils.R')

#' CellMix Command Line Interface
#' 
#' Provides a way to call functions from the \pkg{CellMix} package
#' directly from a console.
#' 
#' @param args command line arguments.
#' This argument is essentially for testing/debugging purposes, since 
#' arguments are usually taken from the command line arguments. 
#' 
CLI <- function(...){

    parser <- CLIArgumentParser()
    
    # commands
    parser$add_command('ged', 'Deconvolve a dataset', default = TRUE)
    parser$add_command('report', 'Generate a deconvolution analysis report')
	
    # parse and run
    parser$parse_cmd()
}

CLI_ged <- function(ARGS = commandArgs(TRUE)){
			
	# define arguments
	parser <- CLIArgumentParser(prog = 'cellmix ged')
	parser$add_argument("dataset", help="Dataset accession number or file to analyse.")
    parser$add_argument("-m", "--method", help="Deconvolution method")
	parser$add_argument("-s", "--sample-groups", help="Name of the variable or file that defines the group of samples to compare.")
	parser$add_argument("-p", "--sample-pairs", help="Name of the variable or file that defines which samples should be paired together.")
	parser$add_argument("-b", "--block-design", help="Name of the variable or file that defines distinct blocks of samples to analyse separately.")
	parser$add_argument("-a", "--annotation-file", help="File that contains variables whose names can be used in other arguments.")
	parser$add_argument("-d", "--data-source", help="Data source where the dataset is fetched from.")
	parser$add_argument("-P", "--proportions", help="Prior proportions to use in the estimation.")
	parser$add_argument("-S", "--signatures", help="Prior cell-specific signatures to use in the estimation.")
	parser$add_argument("-v", "--verbose-level", type = 'integer', default = 1L, help="Verbosity level")
	parser$add_argument("-q", "--quiet", action = 'store_false', dest = 'verbose_level', default = FALSE, help="Run as quietly as possible")
	#
	
	# return parser if no arguments
	if( nargs() == 0L || !length(ARGS) ) return( parser )
	if( is.character(ARGS) ){
		# parse arguments if needed
		ARGS <- parser$parse_args(ARGS)
	}
	
	verbose <- ARGS$verbose_level
    .hasArgument <- .hasArgument(ARGS)
    
	if( verbose > 1 ) message("Calling: ", paste('cellmix', paste0(commandArgs(TRUE), collapse = ' ')))
	if( verbose ) smessage("Loading CellMix ... ")
	library(CellMix, quietly = verbose <= 1)
	if( verbose ) message('OK') 
	
	readFile <- CellMix:::readFile
	lverbose <- CellMix:::lverbose
	vmessage <- CellMix:::vmessage
    NA_Matrix <- CellMix:::NA_Matrix
    smessage <- function(...) vmessage(..., appendLF = FALSE)
	isExpressionSet <- CellMix:::isExpressionSet
	lverbose(verbose)
	
	if( verbose > 1 ){
		vmessage('Arguments:')
		str(ARGS)
	}
	
	# load arguments in current environment
	e <- environment()
	list2env(ARGS, e)
		
	# load data
	if( is.file(dataset) ){ # from file
		smessage("Loading dataset from file '", dataset, "' ... ")
		object <- readFile(dataset)
        
        if( is.matrix(object) ) object <- ExpressionSet(object)
		vmessage('OK')
	}else{ # from data source
		object <- GEDdownload(dataset, verbose = verbose)
	}
    
	# check object
    smessage("Checking dataset ... ")
	if( !isExpressionSet(object) ){
		stop("Invalid dataset object loaded from '", dataset, "': must be an ExpressionSet object [", class(object), ']')
	}
    message('OK')
	
	# load phenotypic file if needed
	if( .hasArgument(annotation_file) ){
		smessage("Loading annotation from file '", annotation_file, "' ... ")
		pheno <- readFile(annotation_file)
		vmessage('OK')
	}
    
    # load known proportions
    if( .hasArgument(proportions) ){
        smessage("Loading proportions from file '", proportions, "' ... ")
        prop <- readFile(proportions)
        vmessage('OK')
        smessage("Checking proportions ... ")
        if( !is.data.frame(prop) ){
            stop('Invalid proportion data: loaded object must be a data.frame [', class(prop), ']')
        }
        # create full proportion matrix
        prop_mat <- NA_Matrix(nrow(prop), ncol(object), dimnames = list(rownames(prop), sampleNames(object)))
        if( is.null(colnames(prop)) ){# number must match
            if( ncol(prop) != ncol(object) ){
                message('ERROR')
                stop("Invalid proportion data: number of columns [", ncol(prop), "] does not match the number of data samples[", ncol(object), "]")
            }
            sn <- 1:ncol(prop)
        }else{
            sn <- intersect(colnames(prop_mat), colnames(prop))
            if( !length(sn) ){
                message('ERROR')
                stop("Invalid proportion data: none of the number of column names match data sample names")
            }
        }
        prop_mat[, sn] <- prop[, sn]
        message('OK [', length(sn), '/', ncol(object), ']')
    }
	
	# 
}


#' CMD report: Deconvolution Report
#' 
#' Generates an HTML deconvolution report for a given dataset.
#' 
#' The report is generated from an Rmd file using the \pkg{knitr} package.
#' 
#' @param ARGS command line arguments
#' 
CLI_report <- function(ARGS = commandArgs(TRUE)){
    
    # define arguments
    parser <- CLIArgumentParser(prog = 'cellmix report')
    parser$add_argument("dataset", help="Dataset accession number or file to analyse.")
    parser$add_argument("-s", "--sample-groups", help="Name of the variable or file that defines the group of samples to compare.")
    parser$add_argument("-p", "--sample-pairs", help="Name of the variable or file that defines which samples should be paired together.")
    parser$add_argument("-b", "--block-design", help="Name of the variable or file that defines distinct blocks of samples to analyse separately.")
    parser$add_argument("-a", "--annotation-file", help="File that contains variables whose names can be used in other arguments.")
    parser$add_argument("-d", "--data-source", help="Data source where the dataset is fetched from.")
    parser$add_argument("-P", "--proportions", help="Prior proportions to use in the estimation.")
    parser$add_argument("-S", "--signatures", help="Prior cell-specific signatures to use in the estimation.")
    parser$add_argument("-v", "--verbose-level", type = 'integer', default = 1L, help="Verbosity level")
    parser$add_argument("-q", "--quiet", action = 'store_false', dest = 'verbose_level', default = FALSE, help="Run as quietly as possible")
    #
    
    # return parser if no arguments
    if( nargs() == 0L || !length(ARGS) ) return( parser )
    if( is.character(ARGS) ){
        # parse arguments if needed
        ARGS <- parser$parse_args(ARGS)
    }
    
    verbose <- ARGS$verbose_level
    .hasArgument <- .hasArgument(ARGS)
    
    if( verbose > 1 ) message("Calling: ", paste('cellmix', paste0(commandArgs(TRUE), collapse = ' ')))
    if( verbose ) smessage("Loading CellMix ... ")
    library(CellMix, quietly = verbose <= 1)
    if( verbose ) message('OK') 
    
    readFile <- CellMix:::readFile
    lverbose <- CellMix:::lverbose
    vmessage <- CellMix:::vmessage
    smessage <- function(...) vmessage(..., appendLF = FALSE)
    isExpressionSet <- CellMix:::isExpressionSet
    lverbose(verbose)
    
    if( verbose > 1 ){
        vmessage('Arguments:')
        str(ARGS)
    }
    
    # load arguments in current environment
    e <- environment()
    list2env(ARGS, e)
    
    # load data
    if( is.file(dataset) ){ # from file
        smessage("Loading dataset from file '", dataset, "' ... ")
        object <- readFile(dataset)
        vmessage('OK')
    }else{ # from data source
        object <- GEDdownload(dataset, verbose = verbose)
    }
    
    # check object
    smessage("Checking dataset ... ")
    if( !isExpressionSet(object) ){
        stop("Invalid dataset object loaded from '", dataset, "': must be an ExpressionSet object [", class(object), ']')
    }
    message('OK')
    
    # load phenotypic file if needed
    if( .hasArgument(annotation_file) ){
        smessage("Loading annotation from file '", annotation_file, "' ... ")
        pheno <- readFile(annotation_file)
        vmessage('OK')
    }
    
    # 
}
