\name{mapGenomeBuilds}

\alias{mapGenomeBuilds}
\alias{genomeBuilds}
\alias{listOrganisms}

\title{Mapping between UCSC and Ensembl Genome Builds}

\description{
  \code{genomeBuilds} lists the available genomes for a given
  species while \code{mapGenomeBuilds} maps between UCSC and Ensemble
  genome builds.
}

\usage{
genomeBuilds(organism, style = c("UCSC", "Ensembl"))

mapGenomeBuilds(genome, style = c("UCSC", "Ensembl"))

listOrganisms()
}

\arguments{
  \item{organism}{A character vector of common names or organism}
  \item{genome}{A character vector of genomes equivalent to UCSC version
    or Ensembl Assemblies}
  \item{style}{A single value equivalent to "UCSC" or "Ensembl"
    specifying the output genome}
}

\details{
  \code{genomeBuilds} lists the currently available genomes for a given list of
  organisms. The genomes can be shown as "UCSC" or "Ensembl" IDs determined
  by \code{style}. \code{organism} must be specified as a character
  vector and match common names (i.e "Dog", "Mouse") or organism name
  (i.e "Homo sapiens", "Mus musculus") . A list of
  available organisms can be shown using \code{listOrganisms()}.

  \code{mapGenomeBuilds} provides a mapping between "UCSC" builds and
  "Ensembl" builds. \code{genome} must be specified as a character
  vector and match either a "UCSC" ID or an "Ensembl"
  Id. \code{genomeBuilds} can be used to get a list of available build
  Ids for a given organism. NA's may be present in the output. This would
  occur when the current genome build removed a previously defined
  genome for an organism.

  In both functions, if \code{style} is not specified, "UCSC" is used as
  default.
}

\value{A data.frame of builds for a given organism or genome in the specified
  \code{style}.  If \code{style == "UCSC"}, ucscID, ucscDate and
  ensemblID are given. If \code{style == "Ensembl"}, ensemblID,
  ensemblVersion, ensemblDate, and ucscID are given. The opposing ID is
  given so that it is possible to distinguish between many-to-one
  mappings.}

\references{
  UCSC genome builds \url{https://genome.ucsc.edu/FAQ/FAQreleases.html}
  Ensembl genome builds \url{http://useast.ensembl.org/info/website/archives/assembly.html}
}

\author{
  Valerie Obenchain \email{Valerie.Obenchain@roswellpark.org} and Lori
  Shepherd \email{Lori.Shepherd@roswellpark.org}}

\examples{

listOrganisms()

genomeBuilds("mouse")
genomeBuilds(c("Mouse", "dog", "human"), style="Ensembl")

mapGenomeBuilds(c("canFam3", "GRCm38", "mm9"))
mapGenomeBuilds(c("canFam3", "GRCm38", "mm9"), style="Ensembl")

}
