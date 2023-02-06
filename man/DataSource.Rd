\name{DataSource}
\alias{DataSource}
\alias{isDataSource}
\title{Gene Expression Data Sources}
\usage{
  DataSource(x, ...)

  isDataSource(x)
}
\arguments{
  \item{x}{data source key or object}

  \item{...}{extra argument to allow extension and passed
  to}
}
\description{
  \code{DataSource} is an S3 generic function that creates
  prototype-like objects whose \code{get} method fetches
  datasets from a given source/repository.

  The \pkg{CellMix} package provides DataSource objects for
  GEO (via the \pkg{GEOquery} package), ArrayExpress (via
  the \pkg{ArrayExpress} package), as well as the internal
  CellMix registry -- which itself fetches datasets from
  multiple repositories using other \code{DataSource}
  objects.

  \code{isDataSource} tests if an object is a
  \code{DataSource}.
}

