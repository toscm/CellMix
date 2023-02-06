# Unit tests for data sources
# 
# Author: Renaud Gaujoux
# Created: Apr 2, 2013
###############################################################################

test.DataSource <- function(){
	
	# constants
	geo_ds <- 'GEO'
	ae_ds <- 'ArrayExpress'
	ds <- c('GEO', 'ArrayExpress', 'CellMix')
	
	# checks
	checkTrue( setequal(DataSource(), ds), "No argument: list all datasources")
	checkIdentical( DataSource(NULL),  NULL, "Argument datasource=NULL: returns NULL.")
	lapply(ds, function(d){
		DS <- DataSource(d)
		.msg <- function(...) paste0("Argument datasource=", d," - ", ...)
		checkTrue( isDataSource(DS), .msg('returns a DataSource object') )
		checkIdentical( DS$name, d, .msg("returns name='", d, "'"))
		checkIdentical( DataSource(d, 'GSE222')$name, d, .msg("valid key provided: returns correct DataSource type (forced value)"))
		checkIdentical( DataSource(d, 'aaa')$name, d, .msg("invalid key provided: returns correct DataSource type (forced value)"))
	})
	checkException( DataSource('aaa'),  "Argument datasource = invalid datasource: throw an error.")
	# GEO
	checkIdentical( DataSource(key = 'GSE124')$name,  geo_ds, "Argument x=GEO access key: returns GEO datasource.")
	# ArrayExpress
	lapply(c('E-TABM-124', 'E-MTAB-12', 'E-MEXP-123', 'E-GEOD-111'), function(aek){
		checkIdentical( DataSource( key = aek)$name,  ae_ds, paste0("Argument x=ArrayExpress access key '", aek, "': returns ArayExpress datasource."))
	})
	TRUE
}