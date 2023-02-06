\name{fetchData}
\alias{fetchData}
\title{Fetching Data from Data Sources}
\usage{
  fetchData(x, ..., datasource = NULL)
}
\arguments{
  \item{x}{dataset accession number (e.g.
  \code{'GSE12345'})}

  \item{...}{extra arguments passed to the datasource's
  \code{get} method.}

  \item{datasource}{datasource where the data is fetched
  from. If \code{NULL}, then all known data sources are
  tried}
}
\description{
  \code{fetchData} retrieves some data from a given data
  source.
}

