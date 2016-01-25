\name{seqlevelsStyle}

\alias{seqlevelsStyle}
\alias{seqlevelsStyle<-}
\alias{seqlevelsStyle,character-method}
\alias{seqlevelsStyle,ANY-method}
\alias{seqlevelsStyle<-,ANY-method}
\alias{seqlevelsStyle<-,character-method}

\alias{genomeStyles}
\alias{extractSeqlevels}
\alias{extractSeqlevelsByGroup}
\alias{mapSeqlevels}
\alias{seqlevelsInGroup}


\title{
  Conveniently rename the seqlevels of an object according to a given style
}

\description{
  The \code{seqlevelsStyle} getter and setter can be used to get the current
  seqlevels style of an object and to rename its seqlevels according to a given
  style.
}

\usage{
seqlevelsStyle(x)
seqlevelsStyle(x) <- value

## Related low-level utilities:
genomeStyles(species)
extractSeqlevels(species, style)
extractSeqlevelsByGroup(species, style, group)
mapSeqlevels(seqnames, style, best.only=TRUE, drop=TRUE)
seqlevelsInGroup(seqnames, group, species, style)
}

\arguments{
  \item{x}{
    The object from/on which to get/set the seqlevels style. \code{x} must
    have a \code{seqlevels} method or be a \code{character} vector.
  }
  \item{value}{
     A single character string that sets the seqlevels style for \code{x}.
  }
  \item{species}{
     The genus and species of the organism in question
     separated by a single space. Don't forget to capitalize the genus.
  }
  \item{style}{
     a character vector with a single element to specify the style.
  }
  \item{group}{ 
     Group can be 'auto' for autosomes, 'sex' for
     sex chromosomes/allosomes, 'circular' for circular chromosomes. The 
     default is 'all' which returns all the chromosomes.
  }
  \item{best.only}{
     if \code{TRUE} (the default), then only the "best"
     sequence renaming maps (i.e. the rows with less NAs) are returned.
  }
  \item{drop}{
     if \code{TRUE} (the default), then a vector is returned
     instead of a matrix when the matrix has only 1 row.
  }
  \item{seqnames}{ 
     a character vector containing the labels attached to the
     chromosomes in a given genome for a given style. For example : For
     \emph{Homo sapiens}, NCBI style - they are
     "1","2","3",...,"X","Y","MT"
  }
}

\details{
    \code{seqlevelsStyle(x)}, \code{seqlevelsStyle(x) <- value}:
    Get the current seqlevels style of an object, or rename its seqlevels
    according to the supplied style.

    \code{genomeStyles}:
    Different organizations have different naming conventions for how they
    name the biologically defined sequence elements (usually chromosomes)
    for each organism they support.  The Seqnames package contains a
    database that defines these different conventions.
 
    genomeStyles() returns the list of all supported seqname mappings, 
    one per supported organism. Each mapping is represented as a data frame 
    with 1 column     per seqname style and 1 row per chromosome name 
    (not all chromosomes of a given organism necessarily belong to the mapping).
 
    genomeStyles(species) returns a data.frame only for the given organism 
    with all its supported seqname mappings. 

    \code{extractSeqlevels}:
    Returns a character vector of the seqnames for a single style and species.

    \code{extractSeqlevelsByGroup}:
    Returns a character vector of the seqnames for a single style and species 
    by group. Group can be 'auto' for autosomes, 'sex' for sex chromosomes/
    allosomes, 'circular' for circular chromosomes. The default is 'all' which 
    returns all the chromosomes. 

    \code{mapSeqlevels}:
    Returns a matrix with 1 column per supplied sequence name and 1 row
    per sequence renaming map compatible with the specified style.
    If \code{best.only} is \code{TRUE} (the default), only the "best"
    renaming maps (i.e. the rows with less NAs) are returned.
 
    \code{seqlevelsInGroup}:
    It takes a character vector along with a group and optional style and 
    species.If group is not specified , it returns "all" or standard/top level 
    seqnames.
    Returns a character vector of seqnames after subsetting for the group 
    specified by the user. See examples for more details.
}

\value{
  For \code{seqlevelsStyle}: A single string containing the style of the
  seqlevels in \code{x}, or a character vector containing the styles of the
  seqlevels in \code{x} if the current style cannot be determined
  unambiguously. Note that this information is not stored in \code{x} but
  inferred by looking up a seqlevels style database stored inside the
  \pkg{GenomeInfoDb} package.
 
  For \code{extractSeqlevels} , \code{extractSeqlevelsByGroup} and 
  \code{seqlevelsInGroup}: A character vector of seqlevels 
  for given supported species and group.
 
  For \code{mapSeqlevels}: A matrix with 1 column per supplied sequence 
  name and 1 row per sequence renaming map compatible with the specified style.

  For \code{genomeStyle}: If species is specified returns a data.frame 
  containg the seqlevels style and its mapping for a given organism. If species
  is not specified, a list is returned with one list per species containing
  the seqlevels style with the corresponding mappings.
}

\author{
  Sonali Arora, Martin Morgan, Marc Carlson, H. Pages
}

\examples{
## ---------------------------------------------------------------------
## seqlevelsStyle() getter and setter
## ---------------------------------------------------------------------

## character vectors:
x <- paste0("chr", 1:5)
seqlevelsStyle(x)
seqlevelsStyle(x) <- "NCBI"
x

## GRanges:
library(GenomicRanges)
gr <- GRanges(rep(c("chr2", "chr3", "chrM"), 2), IRanges(1:6, 10))

seqlevelsStyle(gr)
seqlevelsStyle(gr) <- "NCBI"
gr

seqlevelsStyle(gr)
seqlevelsStyle(gr) <- "dbSNP"
gr

seqlevelsStyle(gr)
seqlevelsStyle(gr) <- "UCSC"
gr

## ---------------------------------------------------------------------
## Related low-level utilities
## ---------------------------------------------------------------------

## Genome styles:
names(genomeStyles())
genomeStyles("Homo_sapiens")
"UCSC" \%in\% names(genomeStyles("Homo_sapiens"))

## Extract seqlevels based on species, style and group:
## The 'group' argument can be 'sex', 'auto', 'circular' or 'all'.

## All:
extractSeqlevels(species="Drosophila_melanogaster", style="Ensembl")

## Sex chromosomes:
extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="sex")

## Autosomes:
extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="auto")


## Identify which seqnames belong to a particular 'group':
newchr <- paste0("chr",c(1:22,"X","Y","M","1_gl000192_random","4_ctg9"))
seqlevelsInGroup(newchr, group="sex")

newchr <- as.character(c(1:22,"X","Y","MT"))
seqlevelsInGroup(newchr, group="all","Homo_sapiens","NCBI")

## Identify which seqnames belong to a species and style:
seqnames <- c("chr1","chr9", "chr2", "chr3", "chr10")
all(seqnames \%in\% extractSeqlevels("Homo_sapiens", "UCSC"))

## Find mapped seqlevelsStyles for exsiting seqnames:
mapSeqlevels(c("chrII", "chrIII", "chrM"), "NCBI")
mapSeqlevels(c("chrII", "chrIII", "chrM"), "Ensembl")

}
