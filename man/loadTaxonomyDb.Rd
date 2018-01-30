\name{loadTaxonomyDb}

\alias{loadTaxonomyDb}
\alias{available.species}

\title{
  Return a data.frame that lists the known taxonomy IDs and their
  corresponding organisms.
}

\description{
  NCBI maintains a collection of unique taxonomy IDs and pairs these
  with associated genus and species designations. This function returns
  the set of pre-processed values that we use to check that something is
  a valid Taxonomy ID (or organism).
}

\usage{
  loadTaxonomyDb()
}

\value{
  A data frame with 1 row per genus/species designation and three columns.
  The 1st column is the taxonomy ID. The second columns is the genus and the
  third is the species name.
}

\author{
  Marc Carlson
}

\examples{
## get the data
taxdb <- loadTaxonomyDb()
tail(taxdb)
## which can then be searched etc.
taxdb[grepl('yoelii', taxdb$species), ]
}

\keyword{manip}
