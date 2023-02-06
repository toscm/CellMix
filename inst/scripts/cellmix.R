#!/usr/bin/env Rscript

# Command line interface to CellMix
#
# usage: cellmix --help
#
# Author: Renaud Gaujoux
# Created: May 9, 2013
###############################################################################

# lookup for main CLI function
pkg <- 'CellMix'
master_cli <- system.file('scripts', 'CLI.R', package = pkg)
if( !nzchar(master_cli) ) master_cli <- paste0('~/projects/', pkg,'/pkg/inst/scripts/CLI.R')
source(master_cli, keep.source = TRUE, chdir = TRUE)
if( !exists('CLI', inherits = FALSE) ){
	stop("Could not start command line interface for package '", pkg, "': main entry point function CLI() not found.")
}

# launch the package's CLI
CLI()
