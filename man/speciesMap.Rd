\name{speciesMap}
\alias{speciesMap}
\alias{validTaxIds}
\alias{specData}
\docType{data}

\author{
  Marc Carlson
  }

  \title{
Species Mapper object and/or list of valid taxonomy Ids
}
\description{
A data frame containing a mapping between NCBI Taxonomy IDs and species
OR a list of valid taxomony IDs (pre-subsetted and cleaned up from that same
data.frame).

The reason why is because sometimes we need the entire frame, but other
times we only need the list of integers, (which is much smaller and
loads much more quickly)

The specData object is an even more complex data.frame that returns
genus and species separately.
}
\usage{
  data(speciesMap)
  data(validTaxIds)
}
\format{
  A data frame with 1416382 observations on the following 2 variables.
  \describe{
    \item{\code{taxon}}{a numeric vector}
    \item{\code{species}}{a character vector}
  }
}
\examples{
data(speciesMap)
sapply(speciesMap, class)   #       taxon     species 
                            #    "integer" "character" 
subset(speciesMap, species=="Homo sapiens")$taxon # [1] 9606
}
\keyword{datasets}
