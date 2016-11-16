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
keepSeqlevels(x, value, pruning.mode=c("error", "coarse", "fine", "tidy"))
dropSeqlevels(x, value, pruning.mode=c("error", "coarse", "fine", "tidy"))
renameSeqlevels(x, value)
restoreSeqlevels(x)
standardChromosomes(x, species=NULL)
keepStandardChromosomes(x, species=NULL,
                        pruning.mode=c("error", "coarse", "fine", "tidy"))
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
  \item{pruning.mode}{
    See \code{?seqinfo} for a description of the pruning modes.
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
gr1 <- keepSeqlevels(gr, "chr1", pruning.mode="coarse")
## Drop 'chr1'. Both 'chr2' and 'chr3' are kept.
gr2 <- dropSeqlevels(gr, "chr1", pruning.mode="coarse")

library(Rsamtools)  # for the ex1.bam file
library(GenomicAlignments)  # for readGAlignments()

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
gal <- readGAlignments(fl)
## If 'value' is named, the names are ignored.
seq2 <- keepSeqlevels(gal, c(foo="seq2"), pruning.mode="coarse")
seqlevels(seq2)

##
## List-like objects:
##

grl0 <- GRangesList(A=GRanges("chr2", IRanges(3:2, 5)),
                    B=GRanges(c("chr2", "chrMT"), IRanges(7:6, 15)),
                    C=GRanges(c("chrY", "chrMT"), IRanges(17:16, 25)),
                    D=GRanges())
## See ?seqinfo for a description of the pruning modes.
keepSeqlevels(grl0, "chr2", pruning.mode="coarse")
keepSeqlevels(grl0, "chr2", pruning.mode="fine")
keepSeqlevels(grl0, "chr2", pruning.mode="tidy")

library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
## Pruning mode "coarse" is particularly well suited on a GRangesList
## object that contains exons grouped by transcript:
ex_by_tx <- exonsBy(txdb, by="tx")
seqlevels(ex_by_tx)
ex_by_tx2 <- keepSeqlevels(ex_by_tx, "chr2L", pruning.mode="coarse")
seqlevels(ex_by_tx2)
## Pruning mode "tidy" is particularly well suited on a GRangesList
## object that contains transcripts grouped by gene:
tx_by_gene <- transcriptsBy(txdb, by="gene")
seqlevels(tx_by_gene)
tx_by_gene2 <- keepSeqlevels(tx_by_gene, "chr2L", pruning.mode="tidy")
seqlevels(tx_by_gene2)

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
stopifnot(identical(gal, gal3))

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

##
## GRanges:
##
gr <- GRanges(c(paste0("chr",c(1:3)), "chr1_gl000191_random",
              "chr1_gl000192_random"), IRanges(1:5, width=3))
gr
keepStandardChromosomes(gr, pruning.mode="coarse")

##
## List-like objects:
##

grl <- GRangesList(GRanges("chr1", IRanges(1:2, 5)),
                   GRanges(c("chr1_GL383519v1_alt", "chr1"), IRanges(5:6, 5)))
## Use pruning.mode="coarse" to drop list elements with mixed seqlevels:
keepStandardChromosomes(grl, pruning.mode="coarse")
## Use pruning.mode="tidy" to keep all list elements with ranges on
## standard chromosomes:
keepStandardChromosomes(grl, pruning.mode="tidy")

##
## The set of standard chromosomes should not be affected by the
## particular seqlevel style currently in use:
##

## NCBI
worm <- GRanges(c("I", "II", "foo", "X", "MT"), IRanges(1:5, width=5))
keepStandardChromosomes(worm, pruning.mode="coarse")

## UCSC
seqlevelsStyle(worm) <- "UCSC"
keepStandardChromosomes(worm, pruning.mode="coarse")

## Ensembl
seqlevelsStyle(worm) <- "Ensembl"
keepStandardChromosomes(worm, pruning.mode="coarse")
}

\keyword{methods}
\keyword{utilities}
