\name{NCBI-utils}

\alias{NCBI-utils}
\alias{find_assembly_ftp_dir}


\title{Utility functions to access NCBI resources}

\description{
  Utility functions to access NCBI resources. Not intended to be used
  directly by the end user.
}

\usage{
find_assembly_ftp_dir(assembly_accession, assembly_name=NA)
}

\arguments{
  \item{assembly_accession}{
    A single string containing a GenBank assembly accession (e.g.
    \code{"GCA_000001405.15"}), a RefSeq assembly accession (e.g.
    \code{"GCF_000001405.26"}), or an NCBI assembly name (e.g.
    \code{"GRCh38"}).
  }
  \item{assembly_name}{
    A single string or \code{NA}.
  }
}

\value{
  A length-2 character vector.
}

\examples{
find_assembly_ftp_dir("GCA_000001405.15")
}

\keyword{manip}

