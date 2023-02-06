\name{Grigoryev-markers}
\alias{Grigoryev-markers}
\title{Grigoryev - Cytometry Antigen Markers}
\description{
  Marker gene list for Human tissues created from the
  \code{\link{markersGrigoryev}} dataset.
}
\details{
  Note this list includes markers for both single and
  multiple cell types. Each kind can be selected using
  argument \code{select} when loading it with
  \code{\link{MarkerList}}. Possible values are:

  \describe{ \item{\code{'all'}}{ all kinds, i.e. single
  and multiple;} \item{\code{'single'}}{ markers for single
  cell type only;} \item{\code{'multiple'}}{ markers for
  multiple cell types only.} }
}
\examples{
# all markers in the list
ml <- MarkerList('Grigoryev')
ml
names(ml)

# only the ones for single cell types
ml <- MarkerList('Grigoryev', select = 'single')
names(ml)
}

