\name{seqinfo}

\alias{seqinfo}
\alias{seqinfo<-}
\alias{seqnames}
\alias{seqnames<-}
\alias{seqlevels}
\alias{seqlevels,ANY-method}
\alias{seqlevels<-}
\alias{seqlevels<-,ANY-method}
\alias{sortSeqlevels}
\alias{sortSeqlevels,character-method}
\alias{sortSeqlevels,ANY-method}
\alias{seqlevelsInUse}
\alias{seqlevelsInUse,Vector-method}
\alias{seqlevelsInUse,CompressedList-method}
\alias{seqlevels0}
\alias{seqlengths}
\alias{seqlengths,ANY-method}
\alias{seqlengths<-}
\alias{seqlengths<-,ANY-method}
\alias{isCircular}
\alias{isCircular,ANY-method}
\alias{isCircular<-}
\alias{isCircular<-,ANY-method}
\alias{genome}
\alias{genome,ANY-method}
\alias{genome<-}
\alias{genome<-,ANY-method}


\title{Accessing/modifying sequence information}

\description{
  A set of generic functions for getting/setting/modifying the sequence
  information stored in an object.
}

\usage{
seqinfo(x)
seqinfo(x,
        new2old=NULL,
        pruning.mode=c("error", "coarse", "fine", "tidy")) <- value

seqnames(x)
seqnames(x) <- value

seqlevels(x)
seqlevels(x,
          pruning.mode=c("error", "coarse", "fine", "tidy")) <- value
sortSeqlevels(x, X.is.sexchrom=NA)
seqlevelsInUse(x)
seqlevels0(x)

seqlengths(x)
seqlengths(x) <- value

isCircular(x)
isCircular(x) <- value

genome(x)
genome(x) <- value
}

\arguments{
  \item{x}{
    The object from/on which to get/set the sequence information.
  }
  \item{new2old}{
    The \code{new2old} argument allows the user to rename, drop, add and/or
    reorder the "sequence levels" in \code{x}.

    \code{new2old} can be \code{NULL} or an integer vector with one element
    per row in \link{Seqinfo} object \code{value} (i.e. \code{new2old} and
    \code{value} must have the same length) describing how the "new" sequence
    levels should be mapped to the "old" sequence levels, that is, how the
    rows in \code{value} should be mapped to the rows in \code{seqinfo(x)}.
    The values in \code{new2old} must be >= 1 and <= \code{length(seqinfo(x))}.
    \code{NA}s are allowed and indicate sequence levels that are being added.
    Old sequence levels that are not represented in \code{new2old} will be
    dropped, but this will fail if those levels are in use (e.g. if \code{x}
    is a \link[GenomicRanges]{GRanges} object with ranges defined on those
    sequence levels) unless a pruning mode is specified via the
    \code{pruning.mode} argument (see below).

    If \code{new2old=NULL}, then sequence levels can only be added to the
    existing ones, that is, \code{value} must have at least as many rows
    as \code{seqinfo(x)} (i.e. \code{length(values) >= length(seqinfo(x))})
    and also \code{seqlevels(values)[seq_len(length(seqlevels(x)))]} must be
    identical to \code{seqlevels(x)}.
  }
  \item{pruning.mode}{
    When some of the seqlevels to drop from \code{x} are in use (i.e. have
    ranges on them), the ranges on these sequences need to be removed before
    the seqlevels can be dropped. We call this \emph{pruning}. The
    \code{pruning.mode} argument controls how to \emph{prune} \code{x}.
    Four pruning modes are currently defined: \code{"error"}, \code{"coarse"},
    \code{"fine"}, and \code{"tidy"}. \code{"error"} is the default. In this
    mode, no pruning is done and an error is raised. The other pruning modes
    do the following:
    \itemize{
        \item \code{"coarse"}: Remove the elements in \code{x} where the
              seqlevels to drop are in use. Typically reduces the length of
              \code{x}. Note that if \code{x} is a list-like object (e.g.
              \link[GenomicRanges]{GRangesList},
              \link[GenomicAlignments]{GAlignmentPairs}, or
              \link[GenomicAlignments]{GAlignmentsList}), then any list
              element in \code{x} where at least one of the sequence levels
              to drop is in use is \emph{fully} removed. In other words,
              when \code{pruning.mode="coarse"}, the \code{seqlevels} setter
              will keep or remove \emph{full list elements} and not try to
              change their content. This guarantees that the exact ranges
              (and their order) inside the individual list elements are
              preserved. This can be a desirable property when the list
              elements represent compound features like exons grouped by
              transcript (stored in a \link[GenomicRanges]{GRangesList}
              object as returned by
              \code{\link[GenomicFeatures]{exonsBy}( , by="tx")}),
              or paired-end or fusion reads, etc...

        \item \code{"fine"}: Supported on list-like objects only. Removes
              the ranges that are on the sequences to drop. This removal is
              done within each list element of the original object \code{x}
              and doesn't affect its length or the order of its list elements.
              In other words, the pruned object is guaranteed to be
              \emph{parallel} to the original object.

        \item \code{"tidy"}: Like the \code{"fine"} pruning above but also
              removes the list elements that become empty as the result of
              the pruning. Note that this pruning mode is particularly well
              suited on a \link[GenomicRanges]{GRangesList} object that
              contains transcripts grouped by gene, as returned by
              \code{\link[GenomicFeatures]{transcriptsBy}( , by="gene")}.
              Finally note that, as a convenience, this pruning mode is
              supported on non list-like objects (e.g.
              \link[GenomicRanges]{GRanges} or
              \link[GenomicAlignments]{GAlignments} objects) and, in this
              case, is equivalent to the \code{"coarse"} mode.
    }
    See the "B. DROP SEQLEVELS FROM A LIST-LIKE OBJECT" section in the examples
    below for an extensive illustration of these pruning modes.
  }
  \item{value}{
    Typically a \link{Seqinfo} object for the \code{seqinfo} setter.

    Either a named or unnamed character vector for the \code{seqlevels}
    setter.

    A vector containing the sequence information to store for the other
    setters.
  }
  \item{X.is.sexchrom}{
    A logical indicating whether X refers to the sexual chromosome
    or to chromosome with Roman Numeral X. If \code{NA}, \code{sortSeqlevels}
    does its best to "guess".
  }
}

\details{
  The \link{Seqinfo} class plays a central role for the functions described
  in this man page because:

  \itemize{
  \item All these functions (except \code{seqinfo}, \code{seqlevelsInUse},
        and \code{seqlevels0}) work on a \link{Seqinfo} object.

  \item For classes that implement it, the \code{seqinfo} getter should
        return a \link{Seqinfo} object.

  \item Default \code{seqlevels}, \code{seqlengths}, \code{isCircular},
        and \code{genome} getters and setters are provided.
        By default, \code{seqlevels(x)} does \code{seqlevels(seqinfo(x))},
        \code{seqlengths(x)} does \code{seqlengths(seqinfo(x))},
        \code{isCircular(x)} does \code{isCircular(seqinfo(x))},
        and \code{genome(x)} does \code{genome(seqinfo(x))}.
        So any class with a \code{seqinfo} getter will have all the above
        getters work out-of-the-box. If, in addition, the class defines
        a \code{seqinfo} setter, then all the corresponding setters will
        also work out-of-the-box.

        Examples of containers that have a \code{seqinfo} getter and setter:
        the \link[GenomicRanges]{GRanges} and \link[GenomicRanges]{GRangesList}
        classes in the \pkg{GenomicRanges} package;
        the \link[SummarizedExperiment]{SummarizedExperiment} class in the
        \pkg{SummarizedExperiment} package;
        the \link[GenomicAlignments]{GAlignments},
        \link[GenomicAlignments]{GAlignmentPairs},
        and \link[GenomicAlignments]{GAlignmentsList} classes in the
        \pkg{GenomicAlignments} package;
        the \link[GenomicFeatures]{TxDb} class in the
        \pkg{GenomicFeatures} package;
        the \link[BSgenome]{BSgenome} class in the
        \pkg{BSgenome} package; etc...
      }

      The \pkg{GenomicRanges} package defines \code{seqinfo} and
      \code{seqinfo<-} methods for these low-level data types:
      \code{List} and \code{IntegerRangesList}. Those objects do not
      have the means to formally store sequence information. Thus,
      the wrappers simply store the \code{Seqinfo} object
      within \code{metadata(x)}. Initially, the metadata
      is empty, so there is some effort to generate a reasonable
      default \code{Seqinfo}. The names of any \code{List} are
      taken as the \code{seqnames}, and the \code{universe} of
      \code{IntegerRangesList} is taken as the \code{genome}.
}

\note{
  The full list of methods defined for a given generic can
  be seen with e.g. \code{showMethods("seqinfo")} or
  \code{showMethods("seqnames")} (for the getters),
  and \code{showMethods("seqinfo<-")} or \code{showMethods("seqnames<-")}
  (for the setters aka \emph{replacement methods}).
  Please be aware that this shows only methods defined in packages
  that are currently attached.
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item The \link{seqlevelsStyle} generic getter and setter.

    \item \link{Seqinfo} objects.

    \item \link[GenomicRanges]{GRanges} and \link[GenomicRanges]{GRangesList}
          objects in the \pkg{GenomicRanges} package.

    \item \link[SummarizedExperiment]{SummarizedExperiment} objects in the
          \pkg{SummarizedExperiment} package.

    \item \link[GenomicAlignments]{GAlignments},
          \link[GenomicAlignments]{GAlignmentPairs},
          and \link[GenomicAlignments]{GAlignmentsList} objects in the
          \pkg{GenomicAlignments} package.

    \item \link[GenomicFeatures]{TxDb} objects in the
          \pkg{GenomicFeatures} package.

    \item \link[BSgenome]{BSgenome} objects in the \pkg{BSgenome} package.

    \item \link{seqlevels-wrappers} for convenience wrappers to the
          \code{seqlevels} getter and setter.

    \item \code{\link{rankSeqlevels}}, on which \code{sortSeqlevels} is
          based.
  }
}

\examples{
## Finding overlaps, comparing, and matching operations between objects
## containing genomic ranges require the objects to have the same
## seqlevels, or they fail. So before one can perform these operations,
## it is often necessary to modify the seqlevels in some of the objects
## involved in the operation so that all the objects have the same
## seqlevels. This is typically done with the seqlevels() setter. It can
## rename, drop, add and reorder seqlevels of an object. Examples below
## show how to mofify the seqlevels of GRanges, GRangesList, and TxDb
## objects but the approach is the same for any object that has seqlevels.

## ---------------------------------------------------------------------
## A. MODIFY THE SEQLEVELS OF A GRanges OBJECT
## ---------------------------------------------------------------------
library(GenomicRanges)

gr <- GRanges(rep(c("chr2", "chr3", "chrM"), 2), IRanges(1:6, 10))

## Add new seqlevels:
seqlevels(gr) <- c("chr1", seqlevels(gr), "chr4")
seqlevels(gr)
seqlevelsInUse(gr)

## Reorder existing seqlevels:
seqlevels(gr) <- rev(seqlevels(gr))
seqlevels(gr)

## Drop all unused seqlevels:
seqlevels(gr) <- seqlevelsInUse(gr)

## Drop some seqlevels in use:
seqlevels(gr, pruning.mode="coarse") <- setdiff(seqlevels(gr), "chr3")
gr

## Rename, add, and reorder the seqlevels all at once:
seqlevels(gr) <- c("chr1", chr2="chr2", chrM="Mitochondrion")
seqlevels(gr)

## ---------------------------------------------------------------------
## B. DROP SEQLEVELS FROM A LIST-LIKE OBJECT
## ---------------------------------------------------------------------

grl0 <- GRangesList(A=GRanges("chr2", IRanges(3:2, 5)),
                    B=GRanges(c("chr2", "chrMT"), IRanges(7:6, 15)),
                    C=GRanges(c("chrY", "chrMT"), IRanges(17:16, 25)),
                    D=GRanges())
grl0

grl1 <- grl0
seqlevels(grl1, pruning.mode="coarse") <- c("chr2", "chr5")
grl1  # grl0[[2]] was fully removed! (even if it had a range on chr2)

## If what is desired is to remove the 2nd range in grl0[[2]] only (i.e.
## the chrMT:6-15 range), or, more generally speaking, to remove the
## ranges within each list element that are located on the seqlevels to
## drop, then use pruning.mode="fine" or pruning.mode="tidy":
grl2 <- grl0
seqlevels(grl2, pruning.mode="fine") <- c("chr2", "chr5")
grl2  # grl0[[2]] not removed, but chrMT:6-15 range removed from it

## Like pruning.mode="fine" but also removes grl0[[3]].
grl3 <- grl0
seqlevels(grl3, pruning.mode="tidy") <- c("chr2", "chr5")
grl3

library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
## Pruning mode "coarse" is particularly well suited on a GRangesList
## object that contains exons grouped by transcript:
ex_by_tx <- exonsBy(txdb, by="tx")
seqlevels(ex_by_tx)
seqlevels(ex_by_tx, pruning.mode="coarse") <- "chr2L"
seqlevels(ex_by_tx)
## Pruning mode "tidy" is particularly well suited on a GRangesList
## object that contains transcripts grouped by gene:
tx_by_gene <- transcriptsBy(txdb, by="gene")
seqlevels(tx_by_gene)
seqlevels(tx_by_gene, pruning.mode="tidy") <- "chr2L"
seqlevels(tx_by_gene)

## ---------------------------------------------------------------------
## C. RENAME THE SEQLEVELS OF A TxDb OBJECT
## ---------------------------------------------------------------------

library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
seqlevels(txdb)

seqlevels(txdb) <- sub("chr", "", seqlevels(txdb))
seqlevels(txdb)

seqlevels(txdb) <- paste0("CH", seqlevels(txdb))
seqlevels(txdb)

seqlevels(txdb)[seqlevels(txdb) == "CHM"] <- "M"
seqlevels(txdb)

## Restore original seqlevels:
seqlevels(txdb) <- seqlevels0(txdb)
seqlevels(txdb)

## ---------------------------------------------------------------------
## D. SORT SEQLEVELS IN "NATURAL" ORDER
## ---------------------------------------------------------------------

sortSeqlevels(c("11", "Y", "1", "10", "9", "M", "2"))

seqlevels <- c("chrXI", "chrY", "chrI", "chrX", "chrIX", "chrM", "chrII")
sortSeqlevels(seqlevels)
sortSeqlevels(seqlevels, X.is.sexchrom=TRUE)
sortSeqlevels(seqlevels, X.is.sexchrom=FALSE)

seqlevels <- c("chr2RHet", "chr4", "chrUextra", "chrYHet",
               "chrM", "chrXHet", "chr2LHet", "chrU",
               "chr3L", "chr3R", "chr2R", "chrX")
sortSeqlevels(seqlevels)

gr <- GRanges()
seqlevels(gr) <- seqlevels
sortSeqlevels(gr)

## ---------------------------------------------------------------------
## E. SUBSET OBJECTS BY SEQLEVELS
## ---------------------------------------------------------------------

tx <- transcripts(txdb)
seqlevels(tx)

## Drop 'M', keep all others.
seqlevels(tx, pruning.mode="coarse") <- seqlevels(tx)[seqlevels(tx) != "M"]
seqlevels(tx)

## Drop all except 'ch3L' and 'ch3R'.
seqlevels(tx, pruning.mode="coarse") <- c("ch3L", "ch3R")
seqlevels(tx)

## ---------------------------------------------------------------------
## F. FINDING METHODS
## ---------------------------------------------------------------------

showMethods("seqinfo")
showMethods("seqinfo<-")

showMethods("seqnames")
showMethods("seqnames<-")

showMethods("seqlevels")
showMethods("seqlevels<-")

if (interactive()) {
  library(GenomicRanges)
  ?`GRanges-class`
}
}

\keyword{methods}
