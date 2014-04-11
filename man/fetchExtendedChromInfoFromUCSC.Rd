\name{fetchExtendedChromInfoFromUCSC}

\alias{fetchExtendedChromInfoFromUCSC}

\title{
  Fetching chromosomes info for some of the UCSC genomes
}

\description{
  Fetch the chromosomes info for some of the UCSC genomes.
  Only supports hg38, mm10, dm3, and sacCer3 at the moment.
}

\usage{
fetchExtendedChromInfoFromUCSC(genome,
        goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
}

\arguments{
  \item{genome}{
    A single string specifying the UCSC genome e.g. \code{"sacCer3"}.
  }
  \item{goldenPath_url}{
    A single string specifying the URL to the UCSC goldenPath location.
    This URL is used internally to build the full URL to the 'chromInfo'
    MySQL dump containing chromosomes information for \code{genome}.
    See Details section below.
  }
}

\details{
  Chromosomes information (e.g. names and lengths) for any UCSC genome
  is stored in the UCSC database in the 'chromInfo' table and is normally
  available as a MySQL dump at:
\preformatted{
  goldenPath_url/<genome>/database/chromInfo.txt.gz
}

  \code{fetchExtendedChromInfoFromUCSC} downloads and imports that table
  into a data frame, and keeps only the \code{UCSC_seqlevels} and
  \code{seqlengths} columns (after renaming them). Then it lookups the
  assembly report at NCBI for the corresponding NCBI assembly (e.g. GRCh38
  assembly for hg38), extracts the seqlevels and GenBank accession numbers
  from the report, matches them to each UCSC seqlevels (using some heuristic),
  and adds them to the returned data frame.
}

\value{
  A data frame with 1 row per seqlevel and with columns:
  \itemize{
    \item UCSC_seqlevels: character vector. No NAs.

    \item NCBI_seqlevels: character vector. Can contain NAs (in which case
          \code{fetchExtendedChromInfoFromUCSC} emits a warning).

    \item GenBank_accns: character vector. Can contain NAs (no warnings).

    \item seqlengths: integer vector. No NAs.
  }

  Note that the rows are not sorted in any particular order.
}

\note{
  Only supports the hg38, mm10, dm3, and sacCer3 genomes at the moment.
  More will come...
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item The \code{\link[GenomicRanges]{seqlevels}} getter and setter
          in the \pkg{GenomicRanges} package.

    \item The \code{\link{seqlevelsStyle}} getter and setter.

    \item The \code{\link[BSgenome]{getBSgenome}} utility in the
          \pkg{BSgenome} package for searching the installed BSgenome
          data packages.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. BASIC EXAMPLE
## ---------------------------------------------------------------------
chrominfo <- fetchExtendedChromInfoFromUCSC("sacCer3")
chrominfo

## ---------------------------------------------------------------------
## B. USE fetchExtendedChromInfoFromUCSC() TO PUT UCSC SEQLEVELS ON THE
##    GRCh38 GENOME
## ---------------------------------------------------------------------

## Load the BSgenome.Hsapiens.NCBI.GRCh38 package:
library(BSgenome)
genome <- getBSgenome("GRCh38")
length(seqlevels(genome))
head(seqlevels(genome), n=30)

## Fetch the extended chromosomes info for the hg38 genome:
hg38_chrominfo <- fetchExtendedChromInfoFromUCSC("hg38")
dim(hg38_chrominfo)
head(hg38_chrominfo, n=30)

## Extract the hg38 seqlevels and put the GRCh38 seqlevels on it as the
## names:
hg38_seqlevels <- setNames(hg38_chrominfo$UCSC_seqlevels,
                           hg38_chrominfo$NCBI_seqlevels)

## Set the hg38 seqlevels on 'genome':
seqlevels(genome) <- hg38_seqlevels[seqlevels(genome)]
head(seqlevels(genome), n=30)

## Sanity check:
hg38_seqlengths <- setNames(hg38_chrominfo$seqlengths,
                            hg38_chrominfo$UCSC_seqlevels)
stopifnot(identical(seqlengths(genome),
                    hg38_seqlengths[seqlevels(genome)]))
}

\keyword{manip}
