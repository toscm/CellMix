\name{gedInput}
\alias{anyRequired}
\alias{gedInput}
\alias{gedIO}
\alias{gedOutput}
\alias{isRequired}
\alias{onlyRequired}
\title{Checking Input and Output Data for GED Algorithms}
\usage{
  gedInput(x)

  gedOutput(x)

  gedIO(x = gedAlgorithm())

  isRequired(type, x)

  onlyRequired(type, x)

  anyRequired(x)
}
\arguments{
  \item{x}{an object that defines a GED algorithm, i.e. a
  \code{GEDalgorithm} object, an access key, a registry
  entry, etc..}

  \item{type}{a character string giving the type of input
  to test, i.e. either \code{'Basis'}, \code{'Coef'} or
  \code{'Marker'}.}
}
\description{
  \code{gedInput} tells which input data is required by a
  given, or all GED algorithm(s).

  Similarly \code{gedOutput} tells the type of output
  returned by each \code{ged} algorithm.

  \code{gedIO} combines the information about input data
  and ouput data in a compact way.

  \code{isRequired}, \code{onlyRequired} and
  \code{anyRequired} return a single logical if \code{x} is
  not missing, or a logical vector with one element per
  registered algorithm otherwise.

  \code{onlyRequired} tells if a given input data is the
  only required data for a GED algorithm to run.

  \code{anyRequired} tells if any input data is required at
  all to run a GED algorithm.
}
\details{
  If \code{x} is not missing \code{gedInput} returns a
  logical vector with one element per type of possible
  input: \code{'Basis'}, \code{'Coef'}, \code{'Marker'}. If
  \code{x} is missing, then it returns a matrix of these
  vectors, with one row per algorithm registered in the
  internal registry.
}
\examples{
# which algorithm requires what
gedInput()
# ask for a given algorithm
gedInput('qprog')
# see type of output
gedOutput()
gedOutput('csSAM')
# input/output combined
gedIO()
gedIO('DSection')
}

