### =========================================================================
### Some low-level utilities to fetch data from the UCSC Genome Browser
### -------------------------------------------------------------------------
###


### NOT exported.
fetch_table_from_UCSC <- function(genome, table,
    colnames=NULL, col2classes=NULL,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    url <- paste(goldenPath.url, genome,
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
fetch_chrom_sizes_from_UCSC <- function(genome,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    col2classes <- c(chrom="character", size="integer", fileName="NULL")
    fetch_table_from_UCSC(genome, "chromInfo",
                          colnames=names(col2classes),
                          col2classes=col2classes,
                          goldenPath.url=goldenPath.url)
}

.cached_chrom_info <- new.env(parent=emptyenv())

.get_chrom_info_from_missing_script <- function(genome,
    assembled.molecules.only=FALSE,
    goldenPath.url=getOption("UCSC.goldenPath.url"),
    recache=FALSE)
{
    ans <- .cached_chrom_info[[genome]]
    if (is.null(ans) || recache) {
        ans <- try(suppressWarnings(
                       fetch_chrom_sizes_from_UCSC(genome,
                           goldenPath.url=goldenPath.url)),
                   silent=TRUE)
        if (inherits(ans, "try-error"))
            stop(wmsg("unknown genome: ", genome))
        oo <- orderSeqlevels(ans[ , "chrom"])
        ans <- S4Vectors:::extract_data_frame_rows(ans, oo)
        ans$is_circular <- make_circ_flags_from_circ_seqs(ans[ , "chrom"])
        .cached_chrom_info[[genome]] <- ans
    }
    if (assembled.molecules.only)
        warning(wmsg("'assembled.molecules' was ignored (don't know what ",
                     "the assembled molecules are for genome ", genome, ")"))
    ans
}

.get_chrom_info_from_script <- function(genome, script_path,
    assembled.molecules.only=FALSE,
    goldenPath.url=getOption("UCSC.goldenPath.url"),
    recache=FALSE)
{
    ## Placeholders. Will actually get defined when we source the script.
    GENOME <- NULL               # Expected to be a single string.
    ASSEMBLED_MOLECULES <- NULL  # Expected to be a character vector with no
                                 # NAs, no empty strings, and no duplicates.
    CIRC_SEQS <- NULL            # Expected to be NULL or a subset of
                                 # ASSEMBLED_MOLECULES.
    GET_CHROM_SIZES <- NULL      # Expected to be NULL (i.e. not defined)
                                 # or a function with 1 argument.
    source(script_path, local=TRUE)

    ## Check script sanity.
    if (!identical(GENOME, genome))
        stop(wmsg(script_path, ": script does not seem ",
                  "to be for genome ", genome))
    stopifnot(is.character(ASSEMBLED_MOLECULES),
              !anyNA(ASSEMBLED_MOLECULES),
              all(nzchar(ASSEMBLED_MOLECULES)),
              !anyDuplicated(ASSEMBLED_MOLECULES))
    if (!is.null(CIRC_SEQS))
        stopifnot(is.character(CIRC_SEQS),
                  !anyDuplicated(CIRC_SEQS),
                  all(CIRC_SEQS %in% ASSEMBLED_MOLECULES))

    ans <- .cached_chrom_info[[genome]]
    if (is.null(ans) || recache) {
        if (is.null(GET_CHROM_SIZES)) {
            ans <- fetch_chrom_sizes_from_UCSC(genome,
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
        ans$is_circular <- make_circ_flags_from_circ_seqs(ans[ , "chrom"],
                                                          CIRC_SEQS)
        .cached_chrom_info[[genome]] <- ans
    }
    if (assembled.molecules.only) {
        i <- seq_along(ASSEMBLED_MOLECULES)
        ans <- S4Vectors:::extract_data_frame_rows(ans, i)
    }
    ans
}

### Returns a 3-column data.frame with columns "chrom" (character), "size"
### (integer), and "is_circular" (logical).
get_chrom_info_from_UCSC <- function(genome,
    assembled.molecules.only=FALSE,
    goldenPath.url=getOption("UCSC.goldenPath.url"),
    recache=FALSE)
{
    if (!isSingleString(genome))
        stop(wmsg("'genome' must be a single string"))
    if (!isTRUEorFALSE(assembled.molecules.only))
        stop(wmsg("'assembled.molecules.only' must be TRUE or FALSE"))
    if (!isSingleString(goldenPath.url))
        stop(wmsg("'goldenPath.url' must be a single string"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))

    script_name <- paste0(genome, ".R")
    script_path <- system.file("UCSC_chrom_info", script_name,
                               package="GenomeInfoDb")
    if (!identical(script_path, "")) {
        ans <- .get_chrom_info_from_script(genome, script_path,
                    assembled.molecules.only=assembled.molecules.only,
                    goldenPath.url=goldenPath.url,
                    recache=recache)
    } else {
        ans <- .get_chrom_info_from_missing_script(genome,
                    assembled.molecules.only=assembled.molecules.only,
                    goldenPath.url=goldenPath.url,
                    recache=recache)
    }
    ans
}

