\name{NCBI-utils}

\alias{NCBI-utils}
\alias{find_NCBI_assembly_ftp_dir}


\title{Utility functions to access NCBI resources}

\description{
  Utility functions to access NCBI resources. Not intended to be used
  directly by the end user.
}

\usage{
find_NCBI_assembly_ftp_dir(assembly_accession, assembly_name=NA)
}

\arguments{
  \item{assembly_accession}{
    A single string containing either a GenBank assembly accession
    (e.g. \code{"GCA_000001405.15"}) or a RefSeq assembly accession
    (e.g. \code{"GCF_000001405.26"}).
  }
  \item{assembly_name}{
    A single string or \code{NA}.
  }
}

\value{
  A length-2 character vector:
  \itemize{
    \item The 1st element in the vector is the URL to the FTP dir,
          without the trailing slash.
    \item The the 2nd element in the vector is the prefix used for
          the names of all the files in the FTP dir.
  }
}

\examples{
find_NCBI_assembly_ftp_dir("GCA_000001405.15")
}

\keyword{manip}

