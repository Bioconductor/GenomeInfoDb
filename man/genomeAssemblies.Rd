\name{genomeAssemblies}
\alias{genomeAssemblies}

\title{
Genome assemblies for multiple organism along with Taxonomy Id.
}
\description{
Genome assemblies for multiple organism along with Taxonomy Id.
}
\usage{
genomeAssemblies()
}
\value{
A data.frame is returned which contains four columns -"Common_Name", 
"Scientific_ID", "UCSC_assembly_ID" and "Taxon_ID".
}
\author{
Sonali Arora
}
\examples{
g <- genomeAssemblies()
head(g)
}


