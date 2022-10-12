\name{Seqinfo-class}
\docType{class}

% Class:
\alias{class:Seqinfo}
\alias{Seqinfo-class}
\alias{Seqinfo}

% Accessors:
\alias{length,Seqinfo-method}
\alias{seqnames,Seqinfo-method}
\alias{seqnames<-,Seqinfo-method}
\alias{names,Seqinfo-method}
\alias{names<-,Seqinfo-method}
\alias{seqlevels,Seqinfo-method}
\alias{seqlevels<-,Seqinfo-method}
\alias{seqlengths,Seqinfo-method}
\alias{seqlengths<-,Seqinfo-method}
\alias{isCircular,Seqinfo-method}
\alias{isCircular<-,Seqinfo-method}
\alias{genome,Seqinfo-method}
\alias{genome<-,Seqinfo-method}

% Subsetting:
\alias{[,Seqinfo-method}

% Coercion:
\alias{as.data.frame.Seqinfo}
\alias{as.data.frame,Seqinfo-method}
\alias{coerce,data.frame,Seqinfo-method}
\alias{coerce,DataFrame,Seqinfo-method}

% Display:
\alias{summary.Seqinfo}
\alias{summary,Seqinfo-method}
\alias{show,Seqinfo-method}

% Combine:
\alias{merge.Seqinfo}
\alias{merge,Seqinfo,missing-method}
\alias{merge,missing,Seqinfo-method}
\alias{merge,Seqinfo,NULL-method}
\alias{merge,NULL,Seqinfo-method}
\alias{merge,Seqinfo,Seqinfo-method}
\alias{intersect,Seqinfo,Seqinfo-method}
\alias{update.Seqinfo}
\alias{update,Seqinfo-method}
\alias{checkCompatibleSeqinfo}

\title{Seqinfo objects}

\description{
  A Seqinfo object is used to store basic information about a set of
  genomic sequences, typically chromosomes (but not necessarily).
}

\details{
  A Seqinfo object has one entry per sequence. Each entry contains the
  following information about the sequence:
  \itemize{
    \item The sequence name (a.k.a. the \emph{seqlevel}) e.g. \code{"chr1"}.
    \item The sequence length.
    \item The sequence \emph{circularity flag}. This is a logical
          indicating whether the sequence is circular (\code{TRUE})
          or linear (\code{FALSE}).
    \item Which genome the sequence belongs to e.g. \code{"hg19"}.
  }
  All entries must contain at least the sequence name. The other information
  is optional. In addition, the \emph{seqnames} in a given Seqinfo object must
  be unique, that is, the object is not allowed to have two entries with the
  same sequence name. In other words, the sequence name is used as the
  \emph{primary key} of a Seqinfo object.

  Note that Seqinfo objects are usually not used as standalone objects
  but are instead typically found inside higher level objects like
  \link[GenomicRanges]{GRanges} or \link[GenomicFeatures]{TxDb} objects.
  These higher level objects will generally provide a \code{seqinfo()}
  accessor for getting/setting their Seqinfo component.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{Seqinfo(seqnames, seqlengths=NA, isCircular=NA, genome=NA)}:
      Create a Seqinfo object and populate it with the supplied data.

      One special form of calling the \code{Seqinfo()} constructor is
      to specify only the \code{genome} argument and set it to the name
      of an NCBI assembly (e.g. \code{Seqinfo(genome="GRCh38.p13")})
      or UCSC genome (e.g. \code{Seqinfo(genome="hg38")}), in which
      case the sequence information is fetched from NCBI or UCSC.
      See Examples section below for some examples.
    }
  }
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is a Seqinfo object.

  \describe{
    \item{}{
      \code{length(x)}:
      Return the number of sequences in \code{x}.
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Get/set the names of the sequences in \code{x}.
      Those names must be non-NA, non-empty and unique.
      They are also called the \emph{sequence levels} or the \emph{keys}
      of the Seqinfo object.

      Note that, in general, the end user should not try to alter the
      sequence levels with \code{seqnames(x) <- value}. The recommended way
      to do this is with \code{seqlevels(x) <- value} as described below.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Same as \code{seqnames(x)} and \code{seqnames(x) <- value}.
    }
    \item{}{
      \code{seqlevels(x)}:
      Same as \code{seqnames(x)}.
    }
    \item{}{
      \code{seqlevels(x) <- value}:
      Can be used to rename, drop, add and/or reorder the sequence levels.
      \code{value} must be either a named or unnamed character vector.
      When \code{value} has names, the names only serve the purpose of
      mapping the new sequence levels to the old ones.
      Otherwise (i.e. when \code{value} is unnamed) this mapping is
      implicitly inferred from the following rules:

      (1) If the number of new and old levels are the same, and if the
          positional mapping between the new and old levels shows that
          some or all of the levels are being renamed, and if the levels
          that are being renamed are renamed with levels that didn't exist
          before (i.e. are not present in the old levels), then
          \code{seqlevels(x) <- value} will just rename the sequence levels.
          Note that in that case the result is the same as with
          \code{seqnames(x) <- value} but it's still recommended to use
          \code{seqlevels(x) <- value} as it is safer.

      (2) Otherwise (i.e. if the conditions for (1) are not satisfied)
          \code{seqlevels(x) <- value} will consider that the sequence
          levels are not being renamed and will just perform
          \code{x <- x[value]}.

      See below for some examples.
    }
    \item{}{
      \code{seqlengths(x)}, \code{seqlengths(x) <- value}:
      Get/set the length for each sequence in \code{x}.
    }
    \item{}{
      \code{isCircular(x)}, \code{isCircular(x) <- value}:
      Get/set the circularity flag for each sequence in \code{x}.
    }
    \item{}{
      \code{genome(x)}, \code{genome(x) <- value}:
      Get/set the genome identifier or assembly name for each sequence
      in \code{x}.
    }
  }
}

\section{Subsetting}{
  In the code snippets below, \code{x} is a Seqinfo object.

  \describe{
    \item{}{
      \code{x[i]}:
      A Seqinfo object can be subsetted only by name i.e. \code{i}
      must be a character vector.
      This is a convenient way to drop/add/reorder the entries
      in a Seqinfo object.

      See below for some examples.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a Seqinfo object.

  \describe{
    \item{}{
      \code{as.data.frame(x)}:
      Turns \code{x} into a data frame.
    }
  }
}

\section{Combining Seqinfo objects}{
  Note that we provide no \code{c()} or \code{rbind()} methods for Seqinfo
  objects. Here is why:

  \code{c()} (like \code{rbind()}) is expected to follow an "appending
  semantic", that is, \code{c(x, y)} is expected to form a new object by
  \emph{appending} the entries in \code{y} to the entries in \code{x},
  thus resulting in an object with \code{length(x) + length(y)} entries.
  The problem with such operation is that it won't be very useful in general,
  because it will tend to break the constraint that the seqnames of a
  Seqinfo object must be unique (primary key).

  So instead, a \code{merge()} method is provided, with a more useful
  semantic. \code{merge(x, y)} does the following:
  \itemize{
    \item If an entry in Seqinfo object \code{x} has the same seqname as an
          entry in Seqinfo object \code{y}, then the 2 entries are fusioned
          together to produce a single entry in the result. This fusion only
          happens if the 2 entries contain compatible information.
    \item If 2 entries cannot be fusioned because they contain incompatible
          information (e.g. different seqlengths or different circularity
          flags), then \code{merge(x, y)} fails with an informative error
          of why \code{x} and \code{y} could not be merged.
  }

  We also implement an \code{update()} method for Seqinfo objects.

  See below for the details.

  In the code snippet below, \code{x}, \code{y}, \code{object}, and
  \code{value}, are Seqinfo objects.

  \describe{
    \item{}{
      \code{merge(x, y, ...)}:
      Merge \code{x} and \code{y} into a single Seqinfo object where the
      keys (i.e. the seqnames) are \code{union(seqnames(x), seqnames(y))}.
      If an entry in \code{y} has the same key as an entry in \code{x}, and
      if the two entries contain compatible information (NA values are treated
      as wildcards i.e. they're compatible with anything), then the two entries
      are merged into a single entry in the result.
      If they cannot be merged (because they contain different seqlengths,
      and/or circularity flags, and/or genome identifiers), then an error
      is raised.
      In addition to check for incompatible sequence information,
      \code{merge(x, y)} also compares \code{seqnames(x)} with
      \code{seqnames(y)} and issues a warning if each of them has names
      not in the other. The purpose of these checks is to try to detect
      situations where the user might be combining or comparing objects
      that use different underlying genomes.

      Note that \code{merge()} can take more than two Seqinfo objects,
      in which case the objects are merged from left to right e.g.
      \preformatted{    merge(x1, x2, x3, x4)}
      is equivalent to
      \preformatted{    merge(merge(merge(x1, x2), x3), x4)}
    }
    \item{}{
      \code{intersect(x, y)}: Finds the intersection between
      two \code{Seqinfo} objects by merging them and subsetting for the
      intersection of their sequence names. This makes it easy to avoid
      warnings about each objects not being a subset of the other one
      during overlap operations.
    }
    \item{}{
      \code{update(object, value)}: Update the entries in Seqinfo object
      \code{object} with the corresponding entries in Seqinfo object
      \code{value}. Note that the seqnames in \code{value} must be a subset
      of the seqnames in \code{object}.
    }
  }

  A convenience wrapper, \code{checkCompatibleSeqinfo()}, is provided
  for checking whether 2 objects have compatible Seqinfo components
  or not. \code{checkCompatibleSeqinfo(x, y)} is equivalent to
  \code{merge(seqinfo(x), seqinfo(y))} so will work on any objects
  \code{x} and \code{y} that support \code{\link{seqinfo}()}.

}

\author{H. Pagès}

\seealso{
  \itemize{
    \item The \code{\link{seqinfo}} getter and setter.

    \item The \code{\link{getChromInfoFromNCBI}} and
          \code{\link{getChromInfoFromUCSC}} utility functions
          that are used behind the scene to generate a Seqinfo
          object for a given assembly/genome (see examples below).
  }
}

\examples{
## ---------------------------------------------------------------------
## A. MAKING A Seqinfo OBJECT FOR A GIVEN NCBI ASSEMBLY OR UCSC GENOME
## ---------------------------------------------------------------------

## One special form of calling the 'Seqinfo()' constructor is to specify
## only the 'genome' argument and set it to the name of an NCBI assembly
## or UCSC genome, in which case the sequence information is fetched
## from NCBI or UCSC ('getChromInfoFromNCBI()' or 'getChromInfoFromUCSC()'
## are used behind the scene for this so internet access is required).

if (interactive()) {
  ## NCBI assemblies (see '?registered_NCBI_assemblies' for the list of
  ## NCBI assemblies that are currently supported):
  Seqinfo(genome="GRCh38")
  Seqinfo(genome="GRCh38.p13")
  Seqinfo(genome="Amel_HAv3.1")
  Seqinfo(genome="WBcel235")
  Seqinfo(genome="TAIR10.1")

  ## UCSC genomes (see '?registered_UCSC_genomes' for the list of UCSC
  ## genomes that are currently supported):
  Seqinfo(genome="hg38")
  Seqinfo(genome="mm10")
  Seqinfo(genome="rn6")
  Seqinfo(genome="bosTau9")
  Seqinfo(genome="canFam3")
  Seqinfo(genome="musFur1")
  Seqinfo(genome="galGal6")
  Seqinfo(genome="dm6")
  Seqinfo(genome="ce11")
  Seqinfo(genome="sacCer3")
}

## ---------------------------------------------------------------------
## B. BASIC MANIPULATION OF A Seqinfo OBJECT
## ---------------------------------------------------------------------

## Note that all the arguments (except 'genome') must have the
## same length. 'genome' can be of length 1, whatever the lengths
## of the other arguments are.
x <- Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chrM"),
             seqlengths=c(100, 200, NA, 15),
             isCircular=c(NA, FALSE, FALSE, TRUE),
             genome="sasquatch")
x

## Accessors:
length(x)
seqnames(x)
names(x)
seqlevels(x)
seqlengths(x)
isCircular(x)
genome(x)

## Get a compact summary:
summary(x)

## Subset by names:
x[c("chrY", "chr3", "chr1")]

## Rename, drop, add and/or reorder the sequence levels:
xx <- x
seqlevels(xx) <- sub("chr", "ch", seqlevels(xx))  # rename
xx
seqlevels(xx) <- rev(seqlevels(xx))  # reorder
xx
seqlevels(xx) <- c("ch1", "ch2", "chY")  # drop/add/reorder
xx
seqlevels(xx) <- c(chY="Y", ch1="1", "22")  # rename/reorder/drop/add
xx

## ---------------------------------------------------------------------
## C. COMBINING 2 Seqinfo OBJECTS
## ---------------------------------------------------------------------

y <- Seqinfo(seqnames=c("chr3", "chr4", "chrM"),
             seqlengths=c(300, NA, 15))
y

## ------ merge() ------

## This issues a warning:
merge(x, y)  # the entries for chr3 and chrM contain information merged
             # from the corresponding entries in 'x' and 'y'

## To get rid of the above warning, either use suppressWarnings() or
## set the genome on 'y':
suppressWarnings(merge(x, y))
genome(y) <- genome(x)
merge(x, y)

## Note that, strictly speaking, merging 2 Seqinfo objects is not
## a commutative operation:
merge(y, x)

## More precisely: In general, 'z1 <- merge(x, y)' is not identical
## to 'z2 <- merge(y, x)'. However 'z1' and 'z2' are guaranteed to
## contain the same information but with their entries possibly in
## different order.

## This contradicts what 'x' says about circularity of chr3 and chrM:
yy <- y
isCircular(yy)[c("chr3", "chrM")] <- c(TRUE, FALSE)

## We say that 'x' and 'yy' are incompatible Seqinfo objects.
yy
if (interactive()) {
  merge(x, yy)  # raises an error
}

## Sanity checks:
stopifnot(identical(x, merge(x, Seqinfo())))
stopifnot(identical(x, merge(Seqinfo(), x)))
stopifnot(identical(x, merge(x, x)))

## ------ update() ------

z <- Seqinfo(seqnames=c("chrM", "chr2", "chr3"),
             seqlengths=c(25, NA, 300),
             genome="chupacabra")
z

update(x, z)

if (interactive()) {
  update(z, x)  # not allowed
  update(x, y)  # not allowed
}

## The seqnames in the 2nd argument can always be forced to be a subset
## of the seqnames in the 1st argument with:
update(x, y[intersect(seqnames(x), seqnames(y))])  # replace entries

## Note that the above is not the same as:
merge(x, y)[seqnames(x)]                           # fusion entries

## The former is guaranteed to work, whatever the Seqinfo objects 'x'
## and 'y'. The latter requires 'x' and 'y' to be compatible.

## Sanity checks:
stopifnot(identical(x, update(x, Seqinfo())))
stopifnot(identical(x, update(x, x)))
stopifnot(identical(z, update(x, z)[seqnames(z)]))

## ---------------------------------------------------------------------
## D. checkCompatibleSeqinfo()
## ---------------------------------------------------------------------
## A simple convenience wrapper to check that 2 objects have compatible
## Seqinfo components.

library(GenomicRanges)
gr1 <- GRanges("chr3:15-25", seqinfo=x)
gr2 <- GRanges("chr3:105-115", seqinfo=y)
if (interactive()) {
  checkCompatibleSeqinfo(gr1, gr2)  # raises an error
}
}

\keyword{methods}
\keyword{classes}
