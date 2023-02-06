# Command Line Interface utils
# 
# Author: Renaud Gaujoux
# Created: May 9, 2013
###############################################################################

.silenceF <- function(f, verbose=FALSE){
    
    if( verbose ) f
    else{
        function(...){
            capture.output(suppressPackageStartupMessages(suppressMessages(res <- f(...)))); 
            res
        }
    }
}

qlibrary <- .silenceF(library, verbose = FALSE)

smessage <- function(...) message(..., appendLF = FALSE)

CLIfile <- function(){
    pattern <- "--file=(.*)"
    if( !length(f <- grep(pattern, commandArgs(), value = TRUE)) ) ''
    else basename(gsub(pattern, "\\1", f))
}

# create a new argument parser
CLIArgumentParser <- function(prog = CLIfile(), ..., epilog = '', show.defaults = TRUE){
    library(argparse, quietly = TRUE)
    
    
    .special <- '__@@##@@__'
    epilog <- paste0(.special, epilog)
    p <- ArgumentParser(prog = prog, ..., epilog = epilog)
    
    # change argument formatter if required
    if( show.defaults ){
        i <- grep("argparse\\.ArgumentParser", p$python_code)
        inst <- p$python_code[i]
        p$python_code[i] <- paste0(substr(inst, 1, nchar(inst)-1), ', formatter_class=argparse.ArgumentDefaultsHelpFormatter)')
    }
    
    p <- proto(p)
    
    # add add_command function
    p$command_loc <- .special
    p$prog <- prog
    p$command <- list()
    p$command_idefault <- 0L
    p$command_default <- function(.){
        if( .$command_idefault ) .$command[.$command_idefault]
        else ''
    }
    
    # add a (sub-)command
    p$add_command <- function(., command, help='', ..., default = FALSE){
        # add command argument if necessary
        if( !length(.$command) ){
            .$.super$add_argument('command', help = paste0(.$prog, ' command to run'))
        }
        # store command
        .$command[command] <- help		
        # store command as default
        if( default ) .$command_idefault <- length(.$command)
    }
    #
    
    p$add_argument <- function(., ..., help = ''){
#		help <- gsub("\n\n", "##NL##", help)
        help <- gsub("\n", "", help)
        .$.super$add_argument(..., help = help)
    }
    
    # overload print_usage
    p$print_usage <- function(.){
        .$.super$print_usage()
        if( length(.$command) ){
            cat("\n  Use --help for listing all available commands\n")
        }
    }
    #
    
    # overload print_help to add command descriptions
    p$print_help <- function(.){
        
        # get formatted help
        h <- paste(capture.output(.$.super$print_help()), collapse="\n")
#		# fix new lines if necessary
#		nl <- strsplit(h, "##NL##")[[1]]
#		if( length(nl) > 1L ){
#			indent <- nchar(gsub("^([ ]+).*", "\\1", tail(strsplit(nl[1], "\n")[[1L]], 1L)))
#			i <- 2:length(nl)
#			print(sprintf(paste0("%", indent, 's'), ''))
#			nl[i] <- paste0(sprintf(paste0("%", indent, 's'), ''), nl[i])
#			h <- paste0(nl, collapse="\n")
#		}
        
        cmds <- ''
        if( length(.$command) ){
            # format command help
            lm <- max(nchar(names(.$command)))
            fmt <- paste0("  %-", lm, "s")
            cmds <- strwrap(.$command, indent = 4, exdent = 2 + lm + 4, width = 80, simplify = FALSE)
            cmds <- sapply(cmds, paste, collapse = "\n")
            cmds <- paste0(sprintf(fmt, names(.$command)), cmds)
            cmds <- paste0('Commands:\n', paste(cmds, collapse = "\n"))
        }
        h <- gsub(.$command_loc, cmds, h, fixed = TRUE)
        cat(h, sep="\n")
    }
    #
    
    # add function call_string
    p$call_string <- function(., args = commandArgs(TRUE)){
        paste(.$prog, paste0(args, collapse = ' '))
    }
    
    # commmand parer
    p$parse_cmd <- function(., ...){
        parseCMD(., ...)
    }
    
    p
}

# combine argument parsers
.combineParser <- function(p1, p2){
    
    if( length(i <- grep("^parser\\.add_argument", p2$python_code)) ){
        p1$.that$python_code <- c(p1$python_code, p2$python_code[i])
    }
    p1
}

.hasArgument <- function(ARGS){
    function(x) length(ARGS[[x]]) && nzchar(ARGS[[x]])
}

logMessage <- function(..., appendLF = TRUE, extfile = NULL){
    
    # output to external file as well
    if( !is.null(extfile) ){
        cat(..., if( appendLF ) "\n", sep ='', file = extfile, append = TRUE)
    }
    message(..., appendLF = appendLF)
    
}

parseCMD <- function(parser, ARGS = commandArgs(TRUE)){
    
    # fix quotes to avoid python JSON parsing error
    ARGS <- gsub("'", "\"", ARGS)
    
    library(pkgmaker, quietly = TRUE)
    # define command line arguments
    prog <- parser$prog
    
    # check validity of command
    # shows usage/help in trivial calls
    if( !length(ARGS) ){
        parser$print_usage()
        return( invisible(parser) )
    }else if( !grepl("^-", ARGS[1L]) ){ # first argument is the command
        command <- ARGS[1L]
        if( !command %in% names(parser$command) ){
            stop("unknown ", prog," command '", command, "'\n"
                    , "  Available commands: ", paste0(names(parser$command), collapse = ', ') 
            #, paste(capture.output(parser$print_usage()), collapse = "\n")
            )
        }
    }else if( any(ARGS %in% c('-h', '--help')) ){
        parser$print_help()
        return( invisible(parser) )
    }else{
        
        # default command if any
        if( nzchar(parser$command_default()) )
            ARGS <- c(parser$command_default(), ARGS)
        else{
            stop("Missing command:\n  "
                    , paste(capture.output(parser$print_usage()), collapse = "\n")
                    , "\n  Available command(s): ", str_out(names(parser$command), Inf, quote=FALSE)
                    , call. = FALSE)
        }
    }
    
    # get command-specific parser
    command <- ARGS[1L]
    if( is.null(cmd_fun <- getFunction(cmd_funame <- paste0('CLI_', command), mustFind = FALSE)) ){
        stop("Could not execute ", prog , " command ", command, ": did not find CLI entry point '", cmd_funame, "'")
    }
    cmd_parser <- cmd_fun(ARGS=NULL)
    ARGS <- ARGS[-1L]
    
    if( !length(ARGS) ){
        # show command line
        cmd_parser$print_usage()
        invisible(cmd_parser)
    }else if( any(ARGS %in% c('-h', '--help')) ){
        cmd_parser$print_help()
        return( invisible(cmd_parser) )
    }else{
        
        # parse command arguments
        args <- cmd_parser$parse_args(ARGS)
        
        # log call and parsed arguments		
        if( FALSE ){
            message('Call: ', parser$call_string(ARGS))
            message('Parsed arguments:')
            str(args)
        }
        #
        
        # call command handler
        cmd_fun(ARGS = args)
    }
}
