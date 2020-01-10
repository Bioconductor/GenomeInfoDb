\name{UCSC-utils}

\alias{get_chrom_info_from_UCSC}

\title{
  Get the chromosome sizes and circularity flags of a UCSC genome assembly
}

\description{
  \code{get_chrom_info_from_UCSC} gets the chromosome sizes and
  circularity flags of a given UCSC genome assembly.
}

\usage{
get_chrom_info_from_UCSC(genome,
                         assembled.molecules.only=FALSE,
                         goldenPath.url=getOption("UCSC.goldenPath.url"),
                         recache=FALSE)
}

\arguments{
  \item{genome}{
    A single string specifying the UCSC genome assembly (e.g.
    \code{"panTro6"}).
  }
  \item{assembled.molecules.only}{
    If \code{FALSE} (the default), sizes and circularity flags are returned
    for all the sequences in the assembly (i.e. all the chromosomes, plasmids,
    and scaffolds).
    If \code{TRUE}, only the sizes and circularity flags of the sequences
    that are considered \emph{assembled molecules} are returned. These are
    the chromosomes (including the mitochondrial chromosome) and plasmids
    only. No scaffolds.
  }
  \item{goldenPath.url}{
    A single string specifying the URL to the UCSC goldenPath location
    where the chromosome sizes are expected to be found.
  }
  \item{recache}{
    \code{get_chrom_info_from_UCSC} uses a cache mechanism so the
    chromosome sizes of a given genome only get downloaded once during
    the current R session (note that the caching is done in memory so
    cached chromosome sizes do NOT persist across sessions).
    Setting \code{recache} to \code{TRUE} forces a new download (and
    recaching) of the chromosome sizes for the specified genome.
  }
}

\value{
  A 3-column data frame with columns \code{chrom} (character),
  \code{size} (integer), and \code{is_circular} (logical).
}

\author{
  H. Pagès
}

\seealso{
  \itemize{
    \item \code{\link{fetchExtendedChromInfoFromUCSC}}
          for fetching more chromosome information for UCSC genomes.

    \item \link{Seqinfo} objects.
  }
}

\examples{
## Internet access required!
get_chrom_info_from_UCSC("hg19")
get_chrom_info_from_UCSC("hg19", assembled.molecules.only=TRUE)
get_chrom_info_from_UCSC("panTro6", assembled.molecules.only=TRUE)
}

\keyword{manip}
