\name{getChromInfoFromNCBI}

\alias{getChromInfoFromNCBI}
\alias{registered_NCBI_assemblies}

\title{
  Get chromosome information for an NCBI assembly
}

\description{
  \code{getChromInfoFromNCBI} returns chromosome information
  like sequence names, lengths and circularity flags for a given
  NCBI assembly e.g. for GRCh38, ARS-UCD1.2, R64, etc...

  Note that \code{getChromInfoFromNCBI} behaves slightly differently
  depending on whether the assembly is \emph{registered} in the
  \pkg{GenomeInfoDb} package or not. See below for the details.

  Use \code{registered_NCBI_assemblies} to list all the NCBI assemblies
  currently registered in the \pkg{GenomeInfoDb} package.
}

\usage{
getChromInfoFromNCBI(assembly,
                     assembled.molecules.only=FALSE,
                     assembly.units=NULL,
                     recache=FALSE,
                     as.Seqinfo=FALSE)

registered_NCBI_assemblies()
}

\arguments{
  \item{assembly}{
    A single string specifying the name of an NCBI assembly (e.g.
    \code{"GRCh38"}).
    Alternatively, an assembly accession (GenBank or RefSeq) can be
    supplied (e.g. \code{"GCF_000001405.12"}).
  }
  \item{assembled.molecules.only}{
    If \code{FALSE} (the default) then chromosome information is returned
    for all the sequences in the assembly (unless \code{assembly.units}
    is specified, see below), that is, for all the chromosomes, plasmids,
    and scaffolds.

    If \code{TRUE} then chromosome information is returned only for the
    \emph{assembled molecules}. These are the chromosomes (including the
    mitochondrial chromosome) and plasmids only. No scaffolds.
  }
  \item{assembly.units}{
    If \code{NULL} (the default) then chromosome information is returned
    for all the sequences in the assembly (unless
    \code{assembled.molecules.only} is set to \code{TRUE}, see above),
    that is, for all the chromosomes, plasmids, and scaffolds.

    \code{assembly.units} can be set to a character vector containing
    the names of \emph{Assembly Units} (e.g. \code{"non-nuclear"}) in
    which case chromosome information is returned only for the sequences
    that belong to these Assembly Units.
  }
  \item{recache}{
    \code{getChromInfoFromNCBI} uses a cache mechanism so the
    chromosome information of a given assembly only gets downloaded
    once during the current R session (note that the caching is done
    in memory so cached information does NOT persist across sessions).
    Setting \code{recache} to \code{TRUE} forces a new download (and
    recaching) of the chromosome information for the specified assembly.
  }
  \item{as.Seqinfo}{
    \code{TRUE} or \code{FALSE} (the default). If \code{TRUE} then a
    \link{Seqinfo} object is returned instead of a data frame. Note that
    only the \code{SequenceName}, \code{SequenceLength}, and \code{circular}
    columns of the data frame are used to make the \link{Seqinfo} object.
    All the other columns are ignored (and lost).
  }
}

\details{
  \emph{registered} vs \emph{unregistered} NCBI assemblies:
  \itemize{
    \item All NCBI assemblies can be looked up by assembly accession
          (GenBank or RefSeq) but only \emph{registered} assemblies
          can also be looked up by assembly name.

    \item For \emph{registered} assemblies, the returned
          circularity flags are guaranteed to be accurate.
          For \emph{unregistered} assemblies, a heuristic is
          used to determine the circular sequences.
  }
  Please contact the maintainer of the \pkg{GenomeInfoDb} package to
  request registration of additional assemblies.
}

\value{
  For \code{getChromInfoFromNCBI}: By default, a 10-column data frame
  with columns:
  \enumerate{
    \item \code{SequenceName}: character.
    \item \code{SequenceRole}: factor.
    \item \code{AssignedMolecule}: factor.
    \item \code{GenBankAccn}: character.
    \item \code{Relationship}: factor.
    \item \code{RefSeqAccn}: character.
    \item \code{AssemblyUnit}: factor.
    \item \code{SequenceLength}: integer. Note that this column **can**
          contain NAs! For example this is the case in assembly Amel_HAv3.1
          where the length of sequence MT is missing or in assembly
          Release 5 where the length of sequence Un is missing.
    \item \code{UCSCStyleName}: character.
    \item \code{circular}: logical.
  }

  For \code{registered_NCBI_assemblies}: A data frame summarizing all the
  NCBI assemblies currently \emph{registered} in the \pkg{GenomeInfoDb}
  package.
}

\author{
  H. Pagès
}

\seealso{
  \itemize{
    \item \code{\link{getChromInfoFromUCSC}} for getting chromosome
          information for a UCSC genome.

    \item \code{\link{getChromInfoFromEnsembl}} for getting chromosome
          information for an Ensembl species.

    \item \link{Seqinfo} objects.
  }
}

\examples{
## Internet access required!

getChromInfoFromNCBI("GRCh37")

getChromInfoFromNCBI("GRCh37", as.Seqinfo=TRUE)

getChromInfoFromNCBI("GRCh37", assembled.molecules.only=TRUE)

getChromInfoFromNCBI("TAIR10.1")

getChromInfoFromNCBI("TAIR10.1", assembly.units="non-nuclear")

## List of NCBI assemblies currently registered in the package:
registered_NCBI_assemblies()

## The GRCh38.p13 assembly only adds "patch sequences" to the GRCh38
## assembly:
GRCh38 <- getChromInfoFromNCBI("GRCh38")
table(GRCh38$SequenceRole)
GRCh38.p13 <- getChromInfoFromNCBI("GRCh38.p13")
table(GRCh38.p13$SequenceRole)  # 185 patch sequences (113 fix + 72 novel)

## Sanity checks:
idx <- match(GRCh38$SequenceName, GRCh38.p13$SequenceName)
stopifnot(!anyNA(idx))
tmp1 <- GRCh38.p13[idx, ]
rownames(tmp1) <- NULL
tmp2 <- GRCh38.p13[-idx, ]
stopifnot(
  identical(tmp1[ , -(5:7)], GRCh38[ , -(5:7)]),
  identical(tmp2, GRCh38.p13[GRCh38.p13$AssemblyUnit == "PATCHES", ])
)
}

\keyword{manip}
