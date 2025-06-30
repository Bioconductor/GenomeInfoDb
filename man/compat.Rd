\name{compat}

\alias{Seqinfo}
\alias{orderSeqlevels}
\alias{rankSeqlevels}

\alias{GenomeDescription}

\title{Functions moved to the Seqinfo package}

\description{
  These functions used to be defined here in the \pkg{GenomeInfoDb} package
  but have moved to the new \pkg{Seqinfo} package in BioC 3.22.
  They are still exported by the \pkg{GenomeInfoDb} package for backward
  compatibility with existing code.

  WARNING: This is a temporary situation only. Packages that \emph{explicitly}
  import these symbols from \pkg{GenomeInfoDb} must be modified to import them
  from the \pkg{Seqinfo} package instead.

  These symbols are actually documented in the \pkg{Seqinfo} package.
  See:
  \itemize{
    \item \code{Seqinfo::\link[Seqinfo]{Seqinfo}}
    \item \code{Seqinfo::\link[Seqinfo]{orderSeqlevels}}
    \item \code{Seqinfo::\link[Seqinfo]{rankSeqlevels}}
    \item \code{Seqinfo::\link[Seqinfo]{GenomeDescription}}
  }
}

\keyword{internal}
