\name{GenomeDescription-class}
\docType{class}

% Classes
\alias{class:GenomeDescription}
\alias{GenomeDescription-class}
\alias{GenomeDescription}

% Accessor methods:
\alias{organism}
\alias{organism,GenomeDescription-method}
\alias{commonName}
\alias{commonName,GenomeDescription-method}
\alias{species}
\alias{species,GenomeDescription-method}
\alias{provider}
\alias{provider,GenomeDescription-method}
\alias{providerVersion}
\alias{providerVersion,GenomeDescription-method}
\alias{releaseDate}
\alias{releaseDate,GenomeDescription-method}
\alias{bsgenomeName}
\alias{bsgenomeName,GenomeDescription-method}
\alias{seqinfo,GenomeDescription-method}
\alias{seqnames,GenomeDescription-method}

% Standard generic methods:
\alias{show,GenomeDescription-method}


\title{GenomeDescription objects}

\description{
  A GenomeDescription object holds the meta information describing
  a given genome.
}

\details{
  In general the user will not need to manipulate directly a
  GenomeDescription instance but will manipulate instead a higher-level
  object that belongs to a class that extends the GenomeDescription class.
  For example, the top-level object defined in any BSgenome data package
  is a \link[BSgenome]{BSgenome} object and the \link[BSgenome]{BSgenome}
  class contains the GenomeDescription class. Thus a \link[BSgenome]{BSgenome}
  object is also a GenomeDescription object and can therefore be treated as
  such. In other words all the methods described below will work on it.
}

\section{Accessor methods}{
  In the code snippets below, \code{object} or \code{x} is a
  GenomeDescription object.

  \describe{
    \item{}{
      \code{organism(object)}:
      Return the scientific name of the organism of the genome e.g.
      \code{"Homo sapiens"}, \code{"Mus musculus"},
      \code{"Caenorhabditis elegans"}, etc...
    }
    \item{}{
      \code{commonName(object)}:
      Return the common name of the organism of the genome e.g.
      \code{"Human"}, \code{"Mouse"}, \code{"Worm"}, etc...
    }
    \item{}{
      \code{provider(x)}:
      Return the provider of this genome e.g. \code{"UCSC"}, \code{"BDGP"},
      \code{"FlyBase"}, etc...
    }
    \item{}{
      \code{releaseDate(x)}:
      Return the release date of this genome e.g. \code{"Mar. 2006"}.
    }
    \item{}{
      \code{bsgenomeName(x)}:
      Uses the meta information stored in \code{x} to make the name of
      the corresponding BSgenome data package (see the
      \code{\link[BSgenome]{available.genomes}} function in the \pkg{BSgenome}
      package for details about the naming scheme used for those packages).
      Of course there is no guarantee that a package with that name
      actually exists.
    }
    \item{\code{seqinfo(x)}}{
      Gets information about the genome sequences. This information
      is returned in a \link{Seqinfo} object.
      Each part of the information can be retrieved separately with
      \code{seqnames(x)}, \code{seqlengths(x)}, and \code{isCircular(x)},
      respectively, as described below.
    }
    \item{\code{seqnames(x)}}{
      Gets the names of the genome sequences.
      \code{seqnames(x)} is equivalent to \code{seqnames(seqinfo(x))}.
    }
    \item{\code{seqlengths(x)}}{
      Gets the lengths of the genome sequences.
      \code{seqlengths(x)} is equivalent to \code{seqlengths(seqinfo(x))}.
    }
    \item{\code{isCircular(x)}}{
      Returns the circularity flags of the genome sequences.
      \code{isCircular(x)} is equivalent to \code{isCircular(seqinfo(x))}.
    }
  }
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item The \code{\link[BSgenome]{available.genomes}} function and the
          \link[BSgenome]{BSgenome} class in the \pkg{BSgenome} package.

    \item The \link{Seqinfo} class.
  }
}

\examples{
library(BSgenome.Celegans.UCSC.ce2)
class(Celegans)
is(Celegans, "GenomeDescription")
provider(Celegans)
seqinfo(Celegans)

genome_desc <- as(Celegans, "GenomeDescription")
class(genome_desc)
genome_desc
provider(genome_desc)
seqinfo(genome_desc)
bsgenomeName(genome_desc)
}

\keyword{methods}
\keyword{classes}
