.onLoad <- function(libname, pkgname)
{
    UCSC.goldenPath.url <- "http://hgdownload.cse.ucsc.edu/goldenPath"
    options(list(UCSC.goldenPath.url=UCSC.goldenPath.url))
}

.test <- function() BiocGenerics:::testPackage("GenomeInfoDb")

