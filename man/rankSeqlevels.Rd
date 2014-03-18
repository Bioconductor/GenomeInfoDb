\name{rankSeqlevels}

\alias{rankSeqlevels}

\title{Assign sequence IDs to sequence names}

\description{
  \code{rankSeqlevels} assigns a unique ID to each unique sequence name
  in the input vector. The returned IDs span 1:N where N is the number of
  unique sequence names in the input vector.
}

\usage{
rankSeqlevels(seqnames, X.is.sexchrom=NA)
}

\arguments{
  \item{seqnames}{
    A character vector or factor containing sequence names.
  }
  \item{X.is.sexchrom}{
    A logical indicating whether X refers to the sexual chromosome 
    or to chromosome with Roman Numeral X. If \code{NA}, \code{rankSeqlevels}
    does its best to "guess".
  }
}

\value{
  An integer vector of the same length as \code{seqnames}.
  The values in the vector span 1:N where N is the number of unique
  sequence names in the input vector.
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item \code{\link[GenomicRanges]{sortSeqlevels}} for sorting the sequence 
    levels of an object in "natural" order.
  }
}

\examples{
library(BSgenome.Scerevisiae.UCSC.sacCer2)
rankSeqlevels(seqnames(Scerevisiae))
rankSeqlevels(seqnames(Scerevisiae)[c(1:5,5:1)])

newchr <- paste0("chr",c(1:3,6:15,4:5,16:22))
newchr
orderSeqlevels(newchr)
rankSeqlevels(newchr)
}

\keyword{manip}

