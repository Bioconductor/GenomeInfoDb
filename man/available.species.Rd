\name{available.species}
\alias{available.species}

\title{
  Returns a data.frame that lists the available species strings and
  their taxonomy Ids.
}

\description{
  NCBI maintains a collection of unique taxonomy Ids and pairs these
  with associated genus and species designations. This function returns
  the set of pre-processed values that we use to check that something is
  a valid Taxonomy ID (or species name)
}

\usage{
  available.species()
}

\value{
  A data frame with 1 row per species designation and two columns.  The
  1st column is the taxonomy Id. The second columns is the species name.
}

\author{
  Marc Carlson
}

\examples{
  ## get the data
  spec <- available.species()
  tail(spec)
  ## which can then be searched etc.
  spec[grepl('yoelii',spec$species),]
}

\keyword{manip}
