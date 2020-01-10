### =========================================================================
### Some low-level utilities to fetch data from the UCSC Genome Browser
### -------------------------------------------------------------------------
###


### NOT exported.
fetch_table_from_UCSC <- function(assembly, table,
    colnames=NULL, col2classes=NULL,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    url <- paste(goldenPath.url, assembly,
                 "database", paste0(table, ".txt.gz"), sep="/")
    destfile <- tempfile()
    download.file(url, destfile, quiet=TRUE)
    if (is.null(col2classes))
        col2classes <- NA
    if (is.null(colnames)) {
        read.table(destfile, sep="\t", quote="",
                             colClasses=col2classes,
                             comment.char="",
                             stringsAsFactors=FALSE)
    } else {
        read.table(destfile, sep="\t", quote="",
                             col.names=colnames,
                             colClasses=col2classes,
                             comment.char="",
                             stringsAsFactors=FALSE)
    }
}

### NOT exported.
fetch_chrom_sizes_from_UCSC <- function(assembly,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    col2classes <- c(chrom="character", size="integer", fileName="NULL")
    fetch_table_from_UCSC(assembly, "chromInfo",
                          colnames=names(col2classes),
                          col2classes=col2classes,
                          goldenPath.url=goldenPath.url)
}

.cached_chrom_sizes <- new.env(parent=emptyenv())

.get_chrom_sizes_from_missing_script <- function(assembly,
    assembled.molecules.only=FALSE,
    goldenPath.url=getOption("UCSC.goldenPath.url"),
    recache=FALSE)
{
    ans <- .cached_chrom_sizes[[assembly]]
    if (is.null(ans) || recache) {
        ans <- try(suppressWarnings(
                       fetch_chrom_sizes_from_UCSC(assembly,
                           goldenPath.url=goldenPath.url)),
                   silent=TRUE)
        if (inherits(ans, "try-error"))
            stop(wmsg("unknown assembly: ", assembly))
        oo <- orderSeqlevels(ans[ , "chrom"])
        ans <- S4Vectors:::extract_data_frame_rows(ans, oo)
        .cached_chrom_sizes[[assembly]] <- ans
    }
    if (assembled.molecules.only)
        warning(wmsg("'assembled.molecules' was ignored (don't ",
                     "know what the assembled molecules are ",
                     "for ", assembly, " assembly)"))
    ans
}

.get_chrom_sizes_from_script <- function(assembly, script_path,
    assembled.molecules.only=FALSE,
    goldenPath.url=getOption("UCSC.goldenPath.url"),
    recache=FALSE)
{
    ## Placeholders. Will actually get defined when we source the script.
    ASSEMBLY <- NULL             # expected to be a single string
    ASSEMBLED_MOLECULES <- NULL  # expected to be a character vector
    GET_CHROM_SIZES <- NULL      # expected to be NULL (i.e. not defined)
                                 # or a function with 1 argument
    source(script_path, local=TRUE)

    ## Check script sanity.
    if (!identical(ASSEMBLY, assembly))
        stop(wmsg(script_path, " script does not seem to be for ",
                  assembly, " assembly"))
    stopifnot(is.character(ASSEMBLED_MOLECULES),
              !anyNA(ASSEMBLED_MOLECULES),
              all(nzchar(ASSEMBLED_MOLECULES)),
              !anyDuplicated(ASSEMBLED_MOLECULES))

    ans <- .cached_chrom_sizes[[assembly]]
    if (is.null(ans) || recache) {
        if (is.null(GET_CHROM_SIZES)) {
            ans <- fetch_chrom_sizes_from_UCSC(assembly,
                                               goldenPath.url=goldenPath.url)
            stopifnot(nrow(ans) == length(ASSEMBLED_MOLECULES))
            oo <- match(ASSEMBLED_MOLECULES, ans[ , "chrom"])
            ans <- S4Vectors:::extract_data_frame_rows(ans, oo)
        } else {
            stopifnot(is.function(GET_CHROM_SIZES))
            ans <- GET_CHROM_SIZES(goldenPath.url=goldenPath.url)
            stopifnot(is.data.frame(ans),
                      identical(sapply(ans, class),
                                c(chrom="character", size="integer")),
                      identical(ans[seq_along(ASSEMBLED_MOLECULES), "chrom"],
                                ASSEMBLED_MOLECULES))
        }
        .cached_chrom_sizes[[assembly]] <- ans
    }
    if (assembled.molecules.only) {
        i <- seq_along(ASSEMBLED_MOLECULES)
        ans <- S4Vectors:::extract_data_frame_rows(ans, i)
    }
    ans
}

### Returns a 2-column data.frame with columns "chrom" (character)
### and "size" (integer).
get_chrom_sizes_from_UCSC <- function(assembly,
    assembled.molecules.only=FALSE,
    goldenPath.url=getOption("UCSC.goldenPath.url"),
    recache=FALSE)
{
    if (!isSingleString(assembly))
        stop(wmsg("'assembly' must be a single string"))
    if (!isTRUEorFALSE(assembled.molecules.only))
        stop(wmsg("'assembled.molecules.only' must be TRUE or FALSE"))
    if (!isSingleString(goldenPath.url))
        stop(wmsg("'goldenPath.url' must be a single string"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))

    script_name <- paste0(assembly, ".R")
    script_path <- system.file("scripts", "UCSC_chrom_sizes", script_name,
                               package="GenomeInfoDb")
    if (!identical(script_path, "")) {
        ans <- .get_chrom_sizes_from_script(assembly, script_path,
                    assembled.molecules.only=assembled.molecules.only,
                    goldenPath.url=goldenPath.url,
                    recache=recache)
    } else {
        ans <- .get_chrom_sizes_from_missing_script(assembly,
                    assembled.molecules.only=assembled.molecules.only,
                    goldenPath.url=goldenPath.url,
                    recache=recache)
    }
    ans
}

