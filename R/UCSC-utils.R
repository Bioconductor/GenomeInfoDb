### =========================================================================
### Some low-level utilities to fetch data from the UCSC Genome Browser
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### UCSC provides a dump of each genome database at:
###   <goldenPath.url>/<genome>/database/
### The database tables are dumped as compressed tab-delimited text files.
### fetch_table_from_UCSC_database() downloads the specified table and
### imports it as a data.frame.
fetch_table_from_UCSC_database <- function(genome, table, col2class,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    url <- paste(goldenPath.url, genome, "database",
                 paste0(table, ".txt.gz"), sep="/")
    fetch_table_from_url(url, colnames=names(col2class), col2class=col2class)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetch_chrom_sizes_from_UCSC()
###

.check_chrom_sizes <- function(chrom_sizes, in_what)
{
    chroms <- chrom_sizes[ , "chrom"]
    if (!is_primary_key(chroms))
        stop(wmsg("invalid data in ", in_what, ": ",
                  "\"chrom\" column contains NAs, empty strings, ",
                  "or duplicates"))
    sizes <- chrom_sizes[ , "size"]
    if (anyNA(sizes) || any(sizes < 0L))
        stop(wmsg("invalid data in ", in_what, ": ",
                  "\"size\" column contains NAs or negative values"))
}

.fetch_chrom_sizes_from_UCSC_database <- function(genome,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    col2class <- c(chrom="character", size="integer", fileName="NULL")
    ans <- fetch_table_from_UCSC_database(genome, "chromInfo",
                                          col2class=col2class,
                                          goldenPath.url=goldenPath.url)
    ## Some sanity checks that should never fail.
    in_what <- paste0("\"chromInfo\" table for UCSC genome ", genome)
    .check_chrom_sizes(ans, in_what)
    ans
}

.fetch_chrom_sizes_from_UCSC_bigZips <- function(genome, use.latest=FALSE,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    filename <- paste0(genome, ".chrom.sizes")
    url <- paste(goldenPath.url, genome, "bigZips", sep="/")
    if (use.latest)
        url <- paste(url, "latest", sep="/")
    url <- paste(url, filename, sep="/")
    col2class <- c(chrom="character", size="integer")
    ans <- fetch_table_from_url(url, colnames=names(col2class),
                                     col2class=col2class)
    ## Some sanity checks that should never fail.
    in_what <- paste0("UCSC file ", filename)
    .check_chrom_sizes(ans, in_what)
    ans
}

### UCSC chrom sizes can be fetched from two locations:
### 1. The 'chromInfo.txt.gz' file at <goldenPath.url>/<genome>/database/:
###    This is a 3-column tab-delimited text file that is compressed.
### 2. The '<genome>.chrom.sizes' file at <goldenPath.url>/<genome>/bigZips/:
###    This is a 2-column tab-delimited text file (uncompressed).
### Except for some rare exceptions, the two files should be available for
### all genomes. Known exceptions (as of Nov 10, 2022):
###   - UCSC does NOT provide 1. for genomes hs1 and mpxvRivers.
###   - UCSC does NOT provide 2. for genomes hg15, rn1, rn2, and rn3.
fetch_chrom_sizes_from_UCSC <- function(genome,
    from=c("database", "bigZips", "bigZips_latest"),
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    if (!isSingleString(genome))
        stop(wmsg("'genome' must be a single string"))
    if (missing(from)) {
        from <- from[[1L]]
    } else {
        if (!isSingleString(from))
            stop(wmsg("'from' must be set to \"database\", ",
                      "\"bigZips\", or \"bigZips_latest\""))
        from <- match.arg(from)
    }
    if (!isSingleString(goldenPath.url))
        stop(wmsg("'goldenPath.url' must be a single string"))

    if (from == "database") {
        ans <- .fetch_chrom_sizes_from_UCSC_database(genome,
                                  goldenPath.url=goldenPath.url)
    } else {
        use_latest <- from == "bigZips_latest"
        ans <- .fetch_chrom_sizes_from_UCSC_bigZips(genome,
                                  use.latest=use_latest,
                                  goldenPath.url=goldenPath.url)
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_UCSC_assembled_molecules_info_table()
###

read_UCSC_assembled_molecules_info_table <- function(file)
{
    col2class <- c(chrom="character", size="integer", circular="logical")
    simple_read_table(file, header=TRUE, col2class=col2class)
}

