\name{seqlevels-wrappers}
\docType{methods}

\alias{seqlevels-wrappers}

\alias{keepSeqlevels}
\alias{dropSeqlevels}
\alias{renameSeqlevels}
\alias{restoreSeqlevels}
\alias{keepStandardChromosomes}

\title{
  Convenience wrappers to the seqlevels() getter and setter
}

\description{
  Keep, drop or rename seqlevels in objects with a \link{Seqinfo} class.
}

\usage{
keepSeqlevels(x, value, ...)
dropSeqlevels(x, value, ...)
renameSeqlevels(x, value, ...)
restoreSeqlevels(x, ...)
keepStandardChromosomes(x, species=NULL)
}

\arguments{
  \item{x}{
    Any object having a \link{Seqinfo} class in which the seqlevels will
    be kept, dropped or renamed.
  }
  \item{value}{
    A named or unnamed character vector.

    Names are ignored by \code{keepSeqlevels} and \code{dropSeqlevels}. 
    Only the values in the character vector dictate which seqlevels to keep
    or drop.

    In the case of \code{renameSeqlevels}, the names are used to map new 
    sequence levels to the old (names correspond to the old levels). When 
    \code{value} is unnamed, the replacement vector must the same length
    and in the same order as the original \code{seqlevels(x)}.
  }
  \item{species}{
    The species name of the \link{Seqinfo} class in which the seqlevels will
    be kept, dropped or renamed.
  }
  \item{...}{
    Arguments passed to other functions.
  }
}

\value{
  The \code{x} object with seqlevels removed or renamed. If \code{x} has
  no seqlevels (empty object) or no replacement values match the current
  seqlevels in \code{x} the unchanged \code{x} is returned.
}

\details{
  Matching and overlap operations on range objects often require that the
  seqlevels match before a comparison can be made (e.g., \code{findOverlaps}).
  \code{keepSeqlevels}, \code{dropSeqlevels} and \code{renameSeqlevels} are 
  high-level convenience functions that wrap the low-level \code{seqlevels}
  function.

  \code{keepSeqlevels}, \code{dropSeqlevels}: Subsetting operations 
  that modify the size of \code{x}. \code{keepSeqlevels} keeps only the
  seqlevels in \code{value} and removes all others. \code{dropSeqlevels}
  drops the levels in \code{value} and retains all others. If \code{value} 
  does not match any seqlevels in \code{x} an empty object is returned.

  \code{renameSeqlevels}: Rename the seqlevels in \code{x} to those in 
  \code{value}. If \code{value} is a named character vector, the names are used
  to map the new seqlevels to the old. When \code{value} is unnamed, the
  replacment vector must be the same length and in the same order as the
  original \code{seqlevels(x)}.

  \code{restoreSeqlevels}: Restore the seqlevels in \code{x} back to the
  original values. Applicable only when \code{x} is a TxDb. The 
  function re-initializes the TxDb which resets the seqlevels, 
  removes masks and any other previous modifications.
  
  \code{keepStandardChromosomes}:Subsetting operation that returns only the 
  'standard' Chromosomes. We define 'standard chromosomes' as those chromosomes 
  which represent sequences in the assembly that are not scaffolds. Also 
  referred to as 'assembly molecule' on NCBI. Applicable when \code{x}
  has a Seqinfo object.This function determines which seqlevels need to be kept 
  using the organism's supported by GenomeInfoDb. The user can also specify the 
  species to get the standard Chromsomes in \code{x}. 
}

\seealso{
  \itemize{
    \item \link{seqinfo} ## Accessing sequence information
    \item \link{Seqinfo} ## The Seqinfo class
  }
}

\examples{

## ---------------------------------------------------------------------
## keepSeqlevels / dropSeqlevels 
## ---------------------------------------------------------------------

## GRanges / GAlignments:

library(GenomicRanges)
gr <- GRanges(c("chr1", "chr1", "chr2", "chr3"), IRanges(1:4, width=3))
seqlevels(gr)
## Keep only 'chr1'
chr1 <- keepSeqlevels(gr, "chr1")
## Drop 'chr1'. Both 'chr2' and 'chr3' are kept.
chr2 <- dropSeqlevels(gr, "chr1")

library(Rsamtools)  # for the ex1.bam file
library(GenomicAlignments)  # for readGAlignments()

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
gal <- readGAlignments(fl)
## If 'value' is named, the names are ignored.
seq2 <- keepSeqlevels(gal, c(foo="seq2"))
seqlevels(seq2)

## GRangesList / GAlignmentsList:

grl <- split(gr, as.character(seqnames(gr)))
dropSeqlevels(grl, c("chr1", "chr2"))
galist <- split(gal, as.character(seqnames(gal)))
keepSeqlevels(galist, "seq2")

## TxDb:

## A TxDb cannot be directly subset with 'keepSeqlevels' 
## and 'dropSeqlevels'.
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
seqlevels(txdb)
\dontrun{
keepSeqlevels(txdb, "chr2L") ## fails
}

## GRanges or GRangesLists extracted from the TxDb can be subset.
txbygene <- transcriptsBy(txdb, "gene")
seqlevels(txbygene)
chr2L <- keepSeqlevels(txbygene, "chr2L")
seqlevels(chr2L)

## ---------------------------------------------------------------------
## renameSeqlevels 
## ---------------------------------------------------------------------

## GAlignments:

seqlevels(gal)
## Rename 'seq2' to 'chr2' with a named or unnamed vector.
gal2a <- renameSeqlevels(gal, c(seq2="chr2"))
gal2b <- renameSeqlevels(gal, c("seq1", "chr2"))
## Names that do not match existing seqlevels are ignored.
## This attempt at renaming does nothing.
gal3 <- renameSeqlevels(gal, c(foo="chr2"))
identical(seqlevels(gal), seqlevels(gal3))

## TxDb:

seqlevels(txdb)
## When the seqlevels of a TxDb are renamed, all future 
## extractions reflect the modified seqlevels.
renameSeqlevels(txdb, sub("chr", "CH", seqlevels(txdb)))
renameSeqlevels(txdb, c(CHM="M"))
seqlevels(txdb)

transcripts <- transcripts(txdb)
identical(seqlevels(txdb), seqlevels(transcripts))

## ---------------------------------------------------------------------
## restoreSeqlevels 
## ---------------------------------------------------------------------

## Restore seqlevels in a TxDb to original values.
\dontrun{
restoreSeqlevels(txdb)
seqlevels(txdb)
}

## ---------------------------------------------------------------------
## keepStandardChromosomes
## ---------------------------------------------------------------------

gr <- GRanges(c(paste0("chr",c(1:3)), "chr1_gl000191_random",
           "chr1_gl000192_random"), IRanges(1:5, width=3))
gr
grl <- split(gr,seqnames(gr))

##GRanges
keepStandardChromosomes(gr)

##GRangesList
keepStandardChromosomes(grl)

plantgr <- GRanges(c(1:5,"MT","Pltd","wrong"), IRanges(1:8,width=5))
keepStandardChromosomes(plantgr, species="Arabidopsis thaliana")
}

\author{Valerie Obenchain \email{vobencha@fhcrc.org}, 
        Sonali Arora \email{sarora@fhcrc.org}}

\keyword{methods}
\keyword{utilities}
