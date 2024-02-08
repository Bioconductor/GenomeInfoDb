.onLoad <- function(libname, pkgname)
{
    UCSC.goldenPath.url <- "https://hgdownload.soe.ucsc.edu/goldenPath"
    options(list(UCSC.goldenPath.url=UCSC.goldenPath.url))
}

.test <- function() BiocGenerics:::testPackage("GenomeInfoDb")

