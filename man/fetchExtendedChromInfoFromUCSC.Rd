\name{fetchExtendedChromInfoFromUCSC}

\alias{fetchExtendedChromInfoFromUCSC}

\title{
  Fetching chromosomes info for some of the UCSC genomes
}

\description{
  Fetch the chromosomes info for some of the UCSC genomes.
  Only supports hg38, hg19, hg18, mm10, mm9, dm3, sacCer3, and sacCer2
  at the moment.
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
  is stored in the UCSC database in the 'chromInfo' table, and is normally
  available as a MySQL dump at:
\preformatted{
  goldenPath_url/<genome>/database/chromInfo.txt.gz
}
  \code{fetchExtendedChromInfoFromUCSC} downloads and imports that table
  into a data frame, keeps only the \code{UCSC_seqlevels} and
  \code{UCSC_seqlengths} columns (after renaming them), and adds the
  \code{circular} logical column.

  Then, if this UCSC genome is based on an NCBI assembly (e.g. hg38 is based
  on GRCh38), the NCBI seqlevels and GenBank accession numbers are extracted
  from the NCBI assembly report and the UCSC seqlevels matched to them (using
  some guided heuristic). Finally the NCBI seqlevels and GenBank accession
  numbers are added to the returned data frame.
}

\value{
  A data frame with 1 row per seqlevel in the UCSC genome, and at least 3
  columns:
  \itemize{
    \item \code{UCSC_seqlevels}: Character vector with no NAs. This is the
          \code{chrom} field of the UCSC 'chromInfo' table for the
          genome. See Details section above.

    \item \code{UCSC_seqlengths}: Integer vector with no NAs. This is the
          \code{size} field of the UCSC 'chromInfo' table for the
          genome. See Details section above.

    \item \code{circular}: Logical vector with no NAs. This knowledge is
          stored in the \pkg{GenomeInfoDb} package itself for the supported
          genomes.
  }
  If the UCSC genome is *not* based on an NCBI assembly (e.g. sacCer2),
  there are no additional columns and a warning is emitted. Note that, in
  this case, the rows are not sorted in any particular order.

  If the UCSC genome is based on an NCBI assembly (e.g. sacCer3),
  the returned data frame has 3 additional columns:
  \itemize{
    \item \code{NCBI_seqlevels}: Character vector. This information is
          obtained from the NCBI assembly report for the genome. Will contain
          NAs for UCSC seqlevels with no corresponding NCBI seqlevels (e.g.
          for chrM in hg18 or chrUextra in dm3), in which case
          \code{fetchExtendedChromInfoFromUCSC} emits a warning.
          
    \item \code{SequenceRole}: Factor with levels \code{assembled-molecule},
          \code{unlocalized-scaffold}, \code{unplaced-scaffold},
          \code{alt-scaffold}, and \code{pseudo-scaffold}. This information
          is obtained from the NCBI assembly report for the genome.
          Can contain NAs but no warning is emitted in that case.

    \item \code{GenBankAccn}: Character vector. This information is obtained
          from the NCBI assembly report for the genome. Can contain NAs but no
          warning is emitted in that case.
  }
  Note that, in this case, the rows are sorted by level in the
  \code{SequenceRole} column, that is, \code{assembled-molecule}s first,
  then \code{unlocalized-scaffold}s, etc, and NAs last.
}

\note{
  Only supports the hg38, hg19, hg18, mm10, mm9, dm3, sacCer3, and sacCer2
  genomes at the moment. More will come...
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

## The sacCer3 UCSC genome is based on an NCBI assembly (RefSeq Assembly
## ID is GCF_000146045.2):
sacCer3_chrominfo <- fetchExtendedChromInfoFromUCSC("sacCer3")
sacCer3_chrominfo

## But the sacCer2 UCSC genome is not:
sacCer2_chrominfo <- fetchExtendedChromInfoFromUCSC("sacCer2")
sacCer2_chrominfo

## ---------------------------------------------------------------------
## B. USING fetchExtendedChromInfoFromUCSC() TO PUT UCSC SEQLEVELS ON
##    THE GRCh38 GENOME
## ---------------------------------------------------------------------

## Load the BSgenome.Hsapiens.NCBI.GRCh38 package:
library(BSgenome)
genome <- getBSgenome("GRCh38")  # this loads the
                                 # BSgenome.Hsapiens.NCBI.GRCh38 package

## A quick look at the GRCh38 seqlevels:
length(seqlevels(genome))
head(seqlevels(genome), n=30)

## Fetch the extended chromosomes info for the hg38 genome:
hg38_chrominfo <- fetchExtendedChromInfoFromUCSC("hg38")
dim(hg38_chrominfo)
head(hg38_chrominfo, n=30)

## 2 sanity checks:
##   1. Check the NCBI seqlevels:
stopifnot(setequal(hg38_chrominfo$NCBI_seqlevels, seqlevels(genome)))
##   2. Check that the sequence lengths in 'hg38_chrominfo' (which are
##      coming from the same 'chromInfo' table as the UCSC seqlevels)
##      are the same as in 'genome':
stopifnot(
  identical(hg38_chrominfo$UCSC_seqlengths,
            unname(seqlengths(genome)[hg38_chrominfo$NCBI_seqlevels]))
)

## Extract the hg38 seqlevels and put the GRCh38 seqlevels on it as
## the names:
hg38_seqlevels <- setNames(hg38_chrominfo$UCSC_seqlevels,
                           hg38_chrominfo$NCBI_seqlevels)

## Set the hg38 seqlevels on 'genome':
seqlevels(genome) <- hg38_seqlevels[seqlevels(genome)]
head(seqlevels(genome), n=30)
}

\keyword{manip}
