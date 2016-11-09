\name{seqlevels-wrappers}
\docType{methods}

\alias{seqlevels-wrappers}

\alias{keepSeqlevels}
\alias{dropSeqlevels}
\alias{renameSeqlevels}
\alias{restoreSeqlevels}
\alias{standardChromosomes}
\alias{keepStandardChromosomes}

\title{
  Convenience wrappers to the seqlevels() getter and setter
}

\description{
  Keep, drop or rename seqlevels in objects with a \link{Seqinfo} class.
}

\usage{
keepSeqlevels(x, value)
dropSeqlevels(x, value)
renameSeqlevels(x, value)
restoreSeqlevels(x)
standardChromosomes(x, species=NULL)
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
     The genus and species of the organism. Supported species can be seen with 
     \code{names(genomeStyles())}.
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
  setter.

  \itemize{
    \item  \code{keepSeqlevels}, \code{dropSeqlevels}: Subsetting operations 
           that modify the size of \code{x}. \code{keepSeqlevels} keeps only
           the seqlevels in \code{value} and removes all others.
           \code{dropSeqlevels} drops the levels in \code{value} and retains
           all others. If \code{value} does not match any seqlevels in \code{x}
           an empty object is returned.

           When \code{x} is a GRangesList it is possible to have 'mixed'
           list elements that have ranges from different chromosomes.
           \code{keepSeqlevels} will not keep 'mixed' list elements 

    \item  \code{renameSeqlevels}: Rename the seqlevels in \code{x} to those in 
           \code{value}. If \code{value} is a named character vector, the names
           are used to map the new seqlevels to the old. When \code{value} is
           unnamed, the replacement vector must be the same length and in the
           same order as the original \code{seqlevels(x)}.

    \item  \code{restoreSeqlevels}: Perform 
           \code{seqlevels(txdb) <- seqlevels0(txdb)}, that is, restore the
           seqlevels in \code{x} back to the original values.
           Applicable only when \code{x} is a TxDb object.

    \item  \code{standardChromosomes}: Lists the 'standard' chromosomes defined
           as sequences in the assembly that are not scaffolds; also referred
           to as an 'assembly molecule' in NCBI.  \code{standardChromosomes}
           attempts to detect the seqlevel style and if more than one style is
           matched, e.g., 'UCSC' and 'Ensembl', the first is chosen.
 
           \code{x} must have a Seqinfo object.  \code{species} can be
           specified as a character string; supported species are listed with
           \code{names(genomeStyles())}.

           When \code{x} contains seqlevels from multiple organisms all
           those considered standard will be kept. For example, if 
           seqlevels are "chr1" and "chr3R" from human and fly both will be
           kept. If \code{species="Homo sapiens"} is specified then only 
           "chr1" is kept.

    \item  \code{keepStandardChromosomes}: Subsetting operation that returns 
           only the 'standard' chromosomes.
 
           \code{x} must have a Seqinfo object.  \code{species} can be
           specified as a character string; supported species are listed with
           \code{names(genomeStyles())}.

           When \code{x} contains seqlevels from multiple organisms all
           those considered standard will be kept. For example, if 
           seqlevels are "chr1" and "chr3R" from human and fly both will be
           kept. If \code{species="Homo sapiens"} is specified then only 
           "chr1" is kept.
    }
}

\seealso{
  \itemize{
    \item \link{seqinfo} ## Accessing sequence information
    \item \link{Seqinfo} ## The Seqinfo class
  }
}

\author{Valerie Obenchain, Sonali Arora}

\examples{

## ---------------------------------------------------------------------
## keepSeqlevels / dropSeqlevels 
## ---------------------------------------------------------------------

##
## GRanges / GAlignments:
##

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

##
## List-like objects
##

grl0 <- GRangesList(GRanges("chr1", IRanges(1:2, 5)),
                    GRanges("chr5", IRanges(11, 18)),
                    GRanges(c("chr4", "chr1"), IRanges(5:6, 5)))
grl0

grl1 <- grl0
keepSeqlevels(grl1, "chr1")
grl1  # grl0[[3]] was fully dropped!

## To keep all ranges within each list element located on "chr2":
grl2 <- grl0[seqnames(grl0) \%in\% "chr1"]
grl2

## The above subsetting doesn't drop any seqlevels:
seqlevels(grl2)

## To drop them:
seqlevels(grl2) <- "chr1"

##
## TxDb:
##

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

##
## GAlignments:
##

seqlevels(gal)
## Rename 'seq2' to 'chr2' with a named vector.
gal2a <- renameSeqlevels(gal, c(seq2="chr2"))
## Rename 'seq2' to 'chr2' with an unnamed vector that includes all 
## seqlevels as they appear in the object.
gal2b <- renameSeqlevels(gal, c("seq1", "chr2"))
## Names that do not match existing seqlevels are ignored.
## This attempt at renaming does nothing.
gal3 <- renameSeqlevels(gal, c(foo="chr2"))
identical(seqlevels(gal), seqlevels(gal3))

##
## TxDb:
##

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
txdb <- restoreSeqlevels(txdb)
seqlevels(txdb)
}

## ---------------------------------------------------------------------
## keepStandardChromosomes
## ---------------------------------------------------------------------

## GRanges
gr <- GRanges(c(paste0("chr",c(1:3)), "chr1_gl000191_random",
              "chr1_gl000192_random"), IRanges(1:5, width=3))
gr
keepStandardChromosomes(gr)

## List-like objects
grl <- GRangesList(GRanges("chr1", IRanges(1:2, 5)),
                   GRanges(c("chr1_GL383519v1_alt", "chr1"), IRanges(5:6, 5)))
## List elements with mixed seqlevels are dropped by default:
keepStandardChromosomes(grl)
## To keep all list elements with ranges on standard chromosomes:
grl[seqnames(grl) \%in\% standardChromosomes(grl)] 

## Keep standard chromosomes of a particular style:
## NCBI
worm <- GRanges(c("I", "II", "foo", "X", "MT"), IRanges(1:5, width=5))
keepStandardChromosomes(worm)

## UCSC
seqlevelsStyle(worm) <- "UCSC"
keepStandardChromosomes(worm)

## Ensembl
seqlevelsStyle(worm) <- "Ensembl"
keepStandardChromosomes(worm)
}

\keyword{methods}
\keyword{utilities}
