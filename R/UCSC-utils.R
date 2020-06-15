### =========================================================================
### Some low-level utilities to fetch data from the UCSC Genome Browser
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


fetch_table_from_UCSC <- function(genome, table, col2class,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    url <- paste(goldenPath.url, genome, "database",
                 paste0(table, ".txt.gz"), sep="/")
    fetch_table_from_url(url, colnames=names(col2class), col2class=col2class)
}

fetch_chrom_sizes_from_UCSC <- function(genome,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    col2class <- c(chrom="character", size="integer", fileName="NULL")
    ans <- fetch_table_from_UCSC(genome, "chromInfo",
                                 col2class=col2class,
                                 goldenPath.url=goldenPath.url)
    ## Should never happen!
    ans_chroms <- ans[ , "chrom"]
    if (!is_primary_key(ans_chroms))
        stop(wmsg("invalid data in \"chromInfo\" table for UCSC genome ",
                  genome, ": \"chrom\" column contains NAs, empty strings, ",
                  "or duplicates"))
    ans_sizes <- ans[ , "size"]
    if (anyNA(ans_sizes) || any(ans_sizes < 0L))
        stop(wmsg("invalid data in \"chromInfo\" table for UCSC genome ",
                  genome, ": \"size\" column contains NAs or negative values"))
    ans
}

