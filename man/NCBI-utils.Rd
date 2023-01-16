\name{NCBI-utils}

\alias{NCBI-utils}
\alias{find_NCBI_assembly_ftp_dir}
\alias{fetch_assembly_report}


\title{Utility functions to access NCBI resources}

\description{
  Low-level utility functions to access NCBI resources. Not intended to be
  used directly by the end user.
}

\usage{
find_NCBI_assembly_ftp_dir(assembly_accession, assembly_name=NA)

fetch_assembly_report(assembly_accession, assembly_name=NA,
                      AssemblyUnits=NULL)
}

\arguments{
  \item{assembly_accession}{
    A single string containing either a GenBank assembly accession
    (e.g. \code{"GCA_000001405.15"}) or a RefSeq assembly accession
    (e.g. \code{"GCF_000001405.26"}).

    Alternatively, for \code{fetch_assembly_report()}, the
    \code{assembly_accession} argument can be set to the URL to
    the \emph{assembly report} (a.k.a. "Full sequence report").
  }
  \item{assembly_name}{
    A single string or \code{NA}.
  }
  \item{AssemblyUnits}{
    By default, all the \emph{assembly units} are included in the data frame
    returned by \code{fetch_assembly_report()}. To include only a subset
    of assembly units, pass a character vector containing the names of the
    assembly units to include to the \code{AssemblyUnits} argument.
  }
}

\value{
  For \code{find_NCBI_assembly_ftp_dir()}: A length-2 character vector:
  \itemize{
    \item The 1st element in the vector is the URL to the FTP dir,
          without the trailing slash.
    \item The 2nd element in the vector is the prefix used in
          the names of most of the files in the FTP dir.
  }

  For \code{fetch_assembly_report()}: A data frame with 1 row per sequence
  in the assembly and 10 columns:
  \enumerate{
    \item SequenceName
    \item SequenceRole
    \item AssignedMolecule
    \item AssignedMoleculeLocationOrType
    \item GenBankAccn
    \item Relationship
    \item RefSeqAccn
    \item AssemblyUnit
    \item SequenceLength
    \item UCSCStyleName
  }
}

\note{
  \code{fetch_assembly_report} is the workhorse behind higher-level
  and more user-friendly \code{\link{getChromInfoFromNCBI}}.
}

\author{
  H. Pagès
}

\seealso{
  \code{\link{getChromInfoFromNCBI}} for a higher-level and
  more user-friendly version of \code{fetch_assembly_report}.
}

\examples{
ftp_dir <- find_NCBI_assembly_ftp_dir("GCA_000001405.15")
ftp_dir

url <- ftp_dir[[1]]     # URL to the FTP dir
prefix <- ftp_dir[[2]]  # prefix used in names of most files

list_ftp_dir(url)

assembly_report_url <- paste0(url, "/", prefix, "_assembly_report.txt")

## To fetch the assembly report for assembly GCA_000001405.15, you can
## call fetch_assembly_report() on the assembly accession or directly
## on the URL to the assembly report:
assembly_report <- fetch_assembly_report("GCA_000001405.15")
dim(assembly_report)
head(assembly_report)

## Sanity check:
assembly_report2 <- fetch_assembly_report(assembly_report_url)
stopifnot(identical(assembly_report, assembly_report2))
}

\keyword{manip}

