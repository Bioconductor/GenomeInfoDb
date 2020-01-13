\name{NCBI-utils}

\alias{get_chrom_info_from_NCBI}
\alias{NCBI_registered_genomes}

\title{
  Get chromosome information of an NCBI genome assembly
}

\description{
  \code{get_chrom_info_from_NCBI} returns chromosome information like
  sequence names, sizes and circularity flags for a given NCBI genome
  assembly e.g. GRCh38, ARS-UCD1.2, R64, etc...

  Note that \code{get_chrom_info_from_NCBI} behaves slightly differently
  depending on whether a genome assembly is \emph{registered} or not in
  the \pkg{GenomeInfoDb} package. See below for the details.

  Use \code{NCBI_registered_genomes} to list all the NCBI genome
  assemblies currently registered in the \pkg{GenomeInfoDb} package.
}

\usage{
get_chrom_info_from_NCBI(genome,
                         assembled.molecules.only=FALSE,
                         recache=FALSE)

NCBI_registered_genomes()
}

\arguments{
  \item{genome}{
    A single string specifying the name of an NCBI genome assembly (e.g.
    \code{"GRCh38"}).
    Alternatively, an assembly accession (GenBank or RefSeq) can be
    supplied (e.g. \code{"GCF_000001405.12"}).
  }
  \item{assembled.molecules.only}{
    If \code{FALSE} (the default) then chromosome information is returned
    for all the sequences in the assembly, that is, for all the chromosomes,
    plasmids, and scaffolds.

    If \code{TRUE} then chromosome information is returned only for the
    \emph{assembled molecules}. These are the chromosomes (including the
    mitochondrial chromosome) and plasmids only. No scaffolds.
  }
  \item{recache}{
    \code{get_chrom_info_from_NCBI} uses a cache mechanism so the
    chromosome information of a given genome assembly only gets downloaded
    once during the current R session (note that the caching is done
    in memory so cached information does NOT persist across sessions).
    Setting \code{recache} to \code{TRUE} forces a new download (and
    recaching) of the chromosome information for the specified genome
    assembly.
  }
}

\details{
  \emph{registered} vs \emph{unregistered} NCBI genome assemblies:
  \itemize{
    \item All NCBI genome assemblies can be looked up by assembly accession
          (GenBank or RefSeq) but only \emph{registered} genome assemblies
          can also be looked up by genome assembly name.

    \item For \emph{registered} genome assemblies, the returned
          circularity flags are guaranteed to be accurate.
          For \emph{unregistered} genome assemblies, a heuristic is
          used to determine the circular sequences.
  }
  Please contact the maintainer of the \pkg{GenomeInfoDb} package to
  request registration of additional genome assemblies.
}

\value{
  \code{get_chrom_info_from_NCBI}: A 9-column data frame with columns:
  \enumerate{
    \item \code{SequenceName} (character).
    \item \code{SequenceRole} (factor).
    \item \code{GenBankAccn} (character).
    \item \code{Relationship} (factor).
    \item \code{RefSeqAccn} (character).
    \item \code{AssemblyUnit} (factor).
    \item \code{SequenceLength} (integer).
    \item \code{UCSCStyleName} (character).
    \item \code{circular} (logical).
  }

  \code{NCBI_registered_genomes}: A data frame summarizing all the NCBI
  genome assemblies currently \emph{registered} in the \pkg{GenomeInfoDb}
  package.
}

\author{
  H. Pagès
}

\seealso{
  \itemize{
    \item \code{\link{get_chrom_info_from_UCSC}} for getting chromosome
          information for a UCSC genome assembly.

    \item \code{\link{fetchExtendedChromInfoFromUCSC}} for fetching extended
          chromosome information of a UCSC genome assembly.

    \item \link{Seqinfo} objects.
  }
}

\examples{
## Internet access required!
get_chrom_info_from_NCBI("GRCh37")
get_chrom_info_from_NCBI("GRCh37", assembled.molecules.only=TRUE)
get_chrom_info_from_NCBI("R64")

## List of NCBI genome assemblies currently registered in the package:
NCBI_registered_genomes()
}

\keyword{manip}
