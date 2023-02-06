\name{hasAnnotation}
\alias{getAnnotation}
\alias{hasAnnotation}
\alias{setAnnotation}
\title{Extracting Annotation from Objects}
\usage{
  hasAnnotation(object, ...)

  getAnnotation(object, ..., null = TRUE)

  setAnnotation(object, value)
}
\arguments{
  \item{object}{an object}

  \item{...}{extra parameters (currently not used)}

  \item{null}{logical that indicates if an empty character
  string should be return as \code{NULL}.}

  \item{value}{new annotation string, e.g.,
  \code{'hgu133plu2.db'}.}
}
\description{
  \code{hasAnnotation} tells if an object has some -- non
  empty -- attached annotation.

  \code{getAnnotation} try extracting embedded annotations.
  By default, it returns \code{NULL} if the object contains
  no annotation.

  \code{setAnnotation} uses a suitable \code{annotation<-}
  method if it exists, or sets the annotation string into
  attribute \code{'annotation'}. It returns the modified
  object.
}

