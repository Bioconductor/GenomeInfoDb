\name{getChromInfoFromUCSC}

\alias{getChromInfoFromUCSC}
\alias{registered_UCSC_genomes}

\title{
  Get chromosome information of a UCSC genome assembly
}

\description{
  \code{getChromInfoFromUCSC} returns chromosome information like
  sequence names, sizes and circularity flags for a given UCSC genome
  assembly e.g. hg19, panTro6, sacCer3, etc...

  Note that \code{getChromInfoFromUCSC} behaves slightly differently
  depending on whether a genome assembly is \emph{registered} or not in
  the \pkg{GenomeInfoDb} package. See below for the details.

  Use \code{registered_UCSC_genomes} to list all the UCSC genome
  assemblies currently registered in the \pkg{GenomeInfoDb} package.
}

\usage{
getChromInfoFromUCSC(genome,
                     assembled.molecules.only=FALSE,
                     goldenPath.url=getOption("UCSC.goldenPath.url"),
                     recache=FALSE)

registered_UCSC_genomes()
}

\arguments{
  \item{genome}{
    A single string specifying the name of a UCSC genome assembly (e.g.
    \code{"panTro6"}).
  }
  \item{assembled.molecules.only}{
    If \code{FALSE} (the default) then chromosome information is returned
    for all the sequences in the assembly, that is, for all the chromosomes,
    plasmids, and scaffolds.

    If \code{TRUE} then chromosome information is returned only for the
    \emph{assembled molecules}. These are the chromosomes (including the
    mitochondrial chromosome) and plasmids only. No scaffolds.

    Note that \code{assembled.molecules.only=TRUE} is supported only for
    \emph{registered} genome assemblies. When used on an \emph{unregistered}
    genome assembly, \code{assembled.molecules.only} is ignored with a warning.
  }
  \item{goldenPath.url}{
    A single string specifying the URL to the UCSC goldenPath location
    where the chromosome sizes are expected to be found.
  }
  \item{recache}{
    \code{getChromInfoFromUCSC} uses a cache mechanism so the
    chromosome sizes of a given genome assembly only get downloaded
    once during the current R session (note that the caching is done
    in memory so cached information does NOT persist across sessions).
    Setting \code{recache} to \code{TRUE} forces a new download (and
    recaching) of the chromosome sizes for the specified genome assembly.
  }
}

\details{
  \emph{registered} vs \emph{unregistered} UCSC genome assemblies:
  \itemize{
    \item \code{assembled.molecules.only=TRUE} is supported only for
          \emph{registered} genome assemblies. For \emph{unregistered}
          genome assemblies, the argument is ignored with a warning.

    \item For \emph{registered} genome assemblies, the returned
          circularity flags are guaranteed to be accurate.
          For \emph{unregistered} genome assemblies, a heuristic is
          used to determine the circular sequences.

    \item For \emph{registered} genome assemblies, special care is taken
          to make sure that the sequences are returned in a sensible order.
          For \emph{unregistered} genome assemblies, a heuristic is used
          to return the sequence in a sensible order.
  }

  Please contact the maintainer of the \pkg{GenomeInfoDb} package to
  request registration of additional genome assemblies.
}

\value{
  \code{getChromInfoFromUCSC}: A 4-column data frame with columns:
  \enumerate{
    \item \code{chrom}: character.
    \item \code{size}: integer.
    \item \code{assembled}: logical.
    \item \code{circular}: logical.
  }

  \code{registered_UCSC_genomes}: A data frame summarizing all the UCSC
  genome assemblies currently \emph{registered} in the \pkg{GenomeInfoDb}
  package.
}

\author{
  H. Pagès
}

\seealso{
  \itemize{
    \item \code{\link{getChromInfoFromNCBI}} for getting chromosome
          information of an NCBI genome assembly.

    \item \code{\link{fetchExtendedChromInfoFromUCSC}} for fetching extended
          chromosome information for a UCSC genome assembly.

    \item \link{Seqinfo} objects.
  }
}

\examples{
## Internet access required!
getChromInfoFromUCSC("hg19")
getChromInfoFromUCSC("hg19", assembled.molecules.only=TRUE)
getChromInfoFromUCSC("panTro6", assembled.molecules.only=TRUE)

## List of UCSC genome assemblies currently registered in the package:
registered_UCSC_genomes()
}

\keyword{manip}
