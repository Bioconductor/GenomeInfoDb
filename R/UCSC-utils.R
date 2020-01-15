### =========================================================================
### Some low-level utilities to fetch data from the UCSC Genome Browser
### -------------------------------------------------------------------------
###


### NOT exported.
fetch_table_from_UCSC <- function(genome, table,
    colnames=NULL, col2class=NULL,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    url <- paste(goldenPath.url, genome,
                 "database", paste0(table, ".txt.gz"), sep="/")
    destfile <- tempfile()
    download.file(url, destfile, quiet=TRUE)
    if (is.null(col2class))
        col2class <- NA
    if (is.null(colnames)) {
        read.table(destfile, sep="\t", quote="",
                             colClasses=col2class,
                             comment.char="",
                             stringsAsFactors=FALSE)
    } else {
        read.table(destfile, sep="\t", quote="",
                             col.names=colnames,
                             colClasses=col2class,
                             comment.char="",
                             stringsAsFactors=FALSE)
    }
}

### NOT exported.
fetch_chrom_sizes_from_UCSC <- function(genome,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    col2class <- c(chrom="character", size="integer", fileName="NULL")
    fetch_table_from_UCSC(genome, "chromInfo",
                          colnames=names(col2class),
                          col2class=col2class,
                          goldenPath.url=goldenPath.url)
}

