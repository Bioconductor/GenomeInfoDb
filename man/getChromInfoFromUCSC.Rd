\name{getChromInfoFromUCSC}

\alias{registered_UCSC_genomes}
\alias{getChromInfoFromUCSC}
\alias{get_and_fix_chrom_info_from_UCSC}
\alias{saveAssembledMoleculesInfoFromUCSC}

\title{
  Get chromosome information for a UCSC genome
}

\description{
  \code{getChromInfoFromUCSC} returns chromosome information like
  sequence names, lengths and circularity flags for a given UCSC genome
  e.g. for hg19, panTro6, sacCer3, etc...

  Note that \code{getChromInfoFromUCSC} behaves slightly differently
  depending on whether a genome is \emph{registered} in the
  \pkg{GenomeInfoDb} package or not. See below for the details.

  Use \code{registered_UCSC_genomes} to list all the UCSC genomes
  currently registered in the \pkg{GenomeInfoDb} package.
}

\usage{
getChromInfoFromUCSC(genome,
                     assembled.molecules.only=FALSE,
                     map.NCBI=FALSE,
                     add.ensembl.col=FALSE,
                     goldenPath.url=getOption("UCSC.goldenPath.url"),
                     recache=FALSE,
                     as.Seqinfo=FALSE)

registered_UCSC_genomes()
}

\arguments{
  \item{genome}{
    A single string specifying the name of a UCSC genome
    (e.g. \code{"panTro6"}).
  }
  \item{assembled.molecules.only}{
    If \code{FALSE} (the default) then chromosome information is returned
    for all the sequences in the genome, that is, for all the chromosomes,
    plasmids, and scaffolds.

    If \code{TRUE} then chromosome information is returned only for the
    \emph{assembled molecules}. These are the chromosomes (including the
    mitochondrial chromosome) and plasmids only. No scaffolds.

    Note that \code{assembled.molecules.only=TRUE} is supported only for
    \emph{registered} genomes. When used on an \emph{unregistered}
    genome, \code{assembled.molecules.only} is ignored with a warning.
  }
  \item{map.NCBI}{
    \code{TRUE} or \code{FALSE} (the default).

    If \code{TRUE} then NCBI chromosome information is bound to
    the result. This information is retrieved from NCBI by calling
    \code{\link{getChromInfoFromNCBI}} on the NCBI assembly that
    the UCSC genome is based on. Then the data frame returned
    by \code{\link{getChromInfoFromNCBI}} ("NCBI chrom info") is
    \emph{mapped} and bound to the data frame returned by
    \code{getChromInfoFromUCSC} ("UCSC chrom info").
    This "map and bind" operation is similar to a \code{JOIN} in SQL.

    Note that not all rows in the "UCSC chrom info" data frame are
    necessarily mapped to a row in the "NCBI chrom info" data frame.
    For example chrM in hg19 has no corresponding sequence in the GRCh37
    assembly (the mitochondrial chromosome was omitted from GRCh37).
    For the unmapped rows the NCBI columns in the final data frame
    are filled with NAs (\code{LEFT JOIN} in SQL).

    The primary use case for using \code{map.NCBI=TRUE} is to map
    UCSC sequence names to NCBI sequence names.
    This is only supported for \emph{registered} UCSC genomes based
    on an NCBI assembly!
  }
  \item{add.ensembl.col}{
    \code{TRUE} or \code{FALSE} (the default). Whether or not the
    Ensembl sequence names should be added to the result (in column
    \code{ensembl}).
  }
  \item{goldenPath.url}{
    A single string specifying the URL to the UCSC goldenPath location
    where the chromosome sizes are expected to be found.
  }
  \item{recache}{
    \code{getChromInfoFromUCSC} uses a cache mechanism so the
    chromosome sizes of a given genome only get downloaded
    once during the current R session (note that the caching is done
    in memory so cached information does NOT persist across sessions).
    Setting \code{recache} to \code{TRUE} forces a new download (and
    recaching) of the chromosome sizes for the specified genome.
  }
  \item{as.Seqinfo}{
    \code{TRUE} or \code{FALSE} (the default). If \code{TRUE} then a
    \link{Seqinfo} object is returned instead of a data frame. Note that
    only the \code{chrom}, \code{size}, and \code{circular} columns of
    the data frame are used to make the \link{Seqinfo} object. All the
    other columns are ignored (and lost).
  }
}

\details{
  *** \emph{Registered} vs \emph{unregistered} UCSC genomes ***
  \itemize{
    \item For \emph{registered} genomes, the returned data frame contains
          information about which sequences are assembled molecules and which
          are not, and the \code{assembled.molecules.only} argument is
          supported. For \emph{unregistered} genomes, this information is
          missing, and the \code{assembled.molecules.only} argument is ignored
          with a warning.

    \item For \emph{registered} genomes, the returned
          circularity flags are guaranteed to be accurate.
          For \emph{unregistered} genomes, a heuristic is
          used to determine the circular sequences.

    \item For \emph{registered} genomes, special care is taken
          to make sure that the sequences are returned in a sensible order.
          For \emph{unregistered} genomes, a heuristic is used
          to return the sequences in a sensible order.
  }

  Please contact the maintainer of the \pkg{GenomeInfoDb} package to
  request registration of additional genomes.

  *** Offline mode ***

  \code{getChromInfoFromUCSC()} supports an "offline mode" when called
  with \code{assembled.molecules.only=TRUE}, but only for a selection of
  \emph{registered} genomes. The "offline mode" works thanks to a collection
  of tab-delimited files stored in the package, that contain the "assembled
  molecules info" for the supported genomes. This makes calls like:
  \preformatted{    getChromInfoFromUCSC("hg38", assembled.molecules.only=TRUE)}
  fast and reliable i.e. the call will always work, even when offline!

  See README.TXT in GenomeInfoDb/inst/extdata/assembled_molecules_info/UCSC/
  for more information.

  Note that calling \code{getChromInfoFromUCSC()} with
  \code{assembled.molecules.only=FALSE} (the default), or with
  \code{recache=TRUE}, will trigger retrieval of the chromosome
  info from UCSC, and will issue a warning if this info no longer
  matches the "assembled molecules info" stored in the package.

  Please contact the maintainer of the \pkg{GenomeInfoDb} package to
  request genome additions to the "offline mode".
}

\value{
  For \code{getChromInfoFromUCSC}: By default, a 4-column data frame
  with columns:
  \enumerate{
    \item \code{chrom}: character.
    \item \code{size}: integer.
    \item \code{assembled}: logical.
    \item \code{circular}: logical.
  }
  If \code{map.NCBI} is \code{TRUE}, then 7 "NCBI columns" are added
  to the result:
  \itemize{
    \item \code{NCBI.SequenceName}: character.
    \item \code{NCBI.SequenceRole}: factor.
    \item \code{NCBI.AssignedMolecule}: factor.
    \item \code{NCBI.GenBankAccn}: character.
    \item \code{NCBI.Relationship}: factor.
    \item \code{NCBI.RefSeqAccn}: character.
    \item \code{NCBI.AssemblyUnit}: factor.
  }
  Note that the names of the "NCBI columns" are those returned
  by \code{\link{getChromInfoFromNCBI}} but with the \code{NCBI.}
  prefix added to them.

  If \code{add.ensembl.col} is \code{TRUE}, the column \code{ensembl}
  is added to the result.

  For \code{registered_UCSC_genomes}: A data frame summarizing all the UCSC
  genomes currently \emph{registered} in the \pkg{GenomeInfoDb} package.
}

\author{
  H. Pagès
}

\seealso{
  \itemize{
    \item \code{\link{getChromInfoFromNCBI}} for getting chromosome
          information for an NCBI assembly.

    \item \code{\link{getChromInfoFromEnsembl}} for getting chromosome
          information for an Ensembl species.

    \item \link{Seqinfo} objects.

    \item The \code{\link[BSgenome]{getBSgenome}} convenience utility in
          the \pkg{BSgenome} package for getting a \link[BSgenome]{BSgenome}
          object from an installed BSgenome data package.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. BASIC EXAMPLES
## ---------------------------------------------------------------------

## --- Internet access required! ---

getChromInfoFromUCSC("hg19")

getChromInfoFromUCSC("hg19", as.Seqinfo=TRUE)

## Map the hg38 sequences to their corresponding sequences in
## the GRCh38.p13 assembly:
getChromInfoFromUCSC("hg38", map.NCBI=TRUE)[c(1, 5)]

## Note that some NCBI-based UCSC genomes contain sequences that
## are not mapped. For example this is the case for chrM in hg19:
hg19 <- getChromInfoFromUCSC("hg19", map.NCBI=TRUE)
hg19[is.na(hg19$NCBI.SequenceName), ]

## Map the hg19 sequences to the Ensembl sequence names:
getChromInfoFromUCSC("hg19", add.ensembl.col=TRUE)

## --- No internet access required! (offline mode) ---

getChromInfoFromUCSC("hg19", assembled.molecules.only=TRUE)

getChromInfoFromUCSC("panTro6", assembled.molecules.only=TRUE)

getChromInfoFromUCSC("bosTau9", assembled.molecules.only=TRUE)

## --- List of UCSC genomes currently registered in the package ---

registered_UCSC_genomes()

## ---------------------------------------------------------------------
## B. USING getChromInfoFromUCSC() TO SET UCSC SEQUENCE NAMES ON THE
##    GRCh38 GENOME
## ---------------------------------------------------------------------

## Load the BSgenome.Hsapiens.NCBI.GRCh38 package:
library(BSgenome)
genome <- getBSgenome("GRCh38")  # this loads the
                                 # BSgenome.Hsapiens.NCBI.GRCh38 package
genome

## Get the chromosome info for the hg38 genome:
hg38_chrom_info <- getChromInfoFromUCSC("hg38", map.NCBI=TRUE)
ncbi2ucsc <- setNames(hg38_chrom_info$chrom,
                      hg38_chrom_info$NCBI.SequenceName)

## Set the UCSC sequence names on 'genome':
seqlevels(genome) <- ncbi2ucsc[seqlevels(genome)]
genome

## Sanity check: check that the sequence lengths in 'genome' are the same
## as in 'hg38_chrom_info':
m <- match(seqlevels(genome), hg38_chrom_info$chrom)
stopifnot(identical(unname(seqlengths(genome)), hg38_chrom_info$size[m]))
}

\keyword{manip}
