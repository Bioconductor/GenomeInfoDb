\name{GenomeInfoDb-helpers}

\alias{GenomeInfoDb-helpers}
\alias{genomeStyles}
\alias{extractSeqlevels}
\alias{extractSeqlevelsByGroup}
\alias{orderSeqlevels}
\alias{mapSeqlevels}
\alias{seqlevelsInGroup}
\alias{seqlevelsStyle}
\alias{seqlevelsStyle,character-method}
\alias{seqlevelsStyle<-}
\alias{seqlevelsStyle<-,character-method}


\title{
List the supported seqname styles for all supported organisms
}
\description{
List the supported seqname styles for all supported organisms
}
\usage{
genomeStyles(species)
extractSeqlevels(species, style)
extractSeqlevelsByGroup(species, style, group)
orderSeqlevels(seqnames, X.is.sexchrom = NA)
mapSeqlevels(seqnames, style, best.only=TRUE, drop=TRUE)
seqlevelsInGroup(seqnames, group, species, style)
seqlevelsStyle(x)
seqlevelsStyle(x) <- value

}
\arguments{
   \item{species}{
      The genus and species of the organism in question
      separated by a single space. Don't forget to capitalize the genus.}
   \item{style}{
      a character vector with a single element to specify the style.}
   \item{group}{ 
      Group can be 'auto' for autosomes, 'sex' for
      sex chromosomes/allosomes, 'linear' for linear chromosomes. The 
      default is 'all' which returns all the chromosomes.}
   \item{best.only}{
      if \code{TRUE} (the default), then only the "best"
      sequence renaming maps (i.e. the rows with less NAs) are returned.}
   \item{drop}{
      if \code{TRUE} (the default), then a vector is returned
      instead of a matrix when the matrix has only 1 row.}
   \item{seqnames}{ 
      a character vector containing the labels attached to the
      chromosomes in a given genome for a given style. For example : For
      \emph{Homo sapiens}, NCBI style - they are
      "1","2","3",...,"X","Y","MT"}
   \item{X.is.sexchrom}{
      A logical indicating whether X refers to the sexual chromosome or
      to chromosome with Roman Numeral X. If \code{NA},
      \code{makeSeqnameIds} does its best to "guess".
   }
   \item{x}{
      The object from/on which to get/set the sequence information.}
   \item{value}{
      A character vector that sets the seqnameStyle for \code{x}.}
}
\value{
  For \code{extractSeqlevels} , \code{extractSeqlevelsByGroup} and 
  \code{seqlevelsInGroup} returns a character vector of seqlevels 
  for given supported species and group.
    
  For \code{mapSeqlevels} returns a matrix with 1 column per supplied sequence 
  name and 1 row per sequence renaming map compatible with the specified style
  
  For \code{seqlevelsStyle} returns a character containing the style of the
  seqlevels supplied. Note that this information is not stored in \code{x}
   but inferred by looking up a seqlevel style database stored inside
   GenomeInfoDb.

  
  For \code{genomeStyle} : If species is specified returns a data.frame 
  containg the seqlevel style and its mapping for a given organism. If species
  is not specified, a list is returned with one list per species containing
  the seqlevel style with the corresponding mappings.
  
  For \code{orderSeqlevels} returns an integer vector with indices of seqlevels
  in their natural order. 
}
\details{
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
    allosomes, 'linear' for linear chromosomes. The default is 'all' which 
    returns all the chromosomes. 
    
    \code{orderSeqlevels}:
    Returns an integer vector while attempting to provide the \dQuote{natural} 
    order of seqnames, e.g.,\code{chr1}, \code{chr2}, \code{chr3}, \ldots
    
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
    
    \code{seqlevelsStyle(x)}: finds the seqlevelsStyle for a given character 
    vector.
}
\author{
    Sonali Arora \email{sarora@fhcrc.org}, Martin Morgan , Marc Carlson, 
    H. Pages
}
\examples{

names(genomeStyles())
genomeStyles("Homo_sapiens")
"UCSC" \%in\% names(genomeStyles("Homo_sapiens"))

## List the supported seqname style for the given species and the given
## style
extractSeqlevels(species="Drosophila_melanogaster" , style="Ensembl")

## List all sex chromosomes for Homo sapiens using style UCSC
## 3 groups are supported: 'auto' for autosomes, 'sex' for allosomes 
## and 'linear' for linear chromosomes  
extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="sex")

## find whether the seqnames belong to a given group
newchr <- paste0("chr",c(1:22,"X","Y","M","1_gl000192_random","4_ctg9"))
seqlevelsInGroup(newchr, group="sex")

newchr <- as.character(c(1:22,"X","Y","MT"))
seqlevelsInGroup(newchr, group="all","Homo_sapiens","NCBI")

## find the seqname Style for a given character vector 
seqlevelsStyle(paste0("chr",c(1:30)))

## order a character vector of seqnames
seqnames <- c("chr1","chr9", "chr2", "chr3", "chr10")
seqnames[orderSeqlevels(seqnames)]

## if we have a vector conatining seqnames and we want to verify the
## species and style for them , we can use:
all(seqnames \%in\% extractSeqlevels("Homo_sapiens", "UCSC"))

## find mapped seqlevelsStyles for exsiting seqnames
mapSeqlevels(c("chrII", "chrIII", "chrM"), "NCBI")
mapSeqlevels(c("chrII", "chrIII", "chrM"), "Ensembl")

}
