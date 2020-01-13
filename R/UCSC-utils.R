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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### UCSC_registered_genomes()
###

.parse_script_for_UCSC_registered_genome <- function(script_path)
{
    ## Placeholders. Will actually get defined when we source the script.
    GENOME <- NULL               # Expected to be a single string.
    ORGANISM <- NULL             # Expected to be a single string.
    ASSEMBLED_MOLECULES <- NULL  # Expected to be a character vector with no
                                 # NAs, no empty strings, and no duplicates.
    CIRC_SEQS <- NULL            # Expected to be NULL or a subset of
                                 # ASSEMBLED_MOLECULES.
    GET_CHROM_SIZES <- NULL      # Expected to be NULL (i.e. not defined)
                                 # or a function with 1 argument.

    source(script_path, local=TRUE)

    ## Check script sanity.

    stopifnot(isSingleString(GENOME))

    stopifnot(isSingleString(ORGANISM))

    ASSEMBLED_MOLECULES <- ASSEMBLED_MOLECULES
    stopifnot(is.character(ASSEMBLED_MOLECULES),
              !anyNA(ASSEMBLED_MOLECULES),
              all(nzchar(ASSEMBLED_MOLECULES)),
              !anyDuplicated(ASSEMBLED_MOLECULES))

    CIRC_SEQS <- CIRC_SEQS
    if (!is.null(CIRC_SEQS))
        stopifnot(is.character(CIRC_SEQS),
                  !anyDuplicated(CIRC_SEQS),
                  all(CIRC_SEQS %in% ASSEMBLED_MOLECULES))

    GET_CHROM_SIZES <- GET_CHROM_SIZES
    if (!is.null(GET_CHROM_SIZES))
        stopifnot(is.function(GET_CHROM_SIZES))

    list(GENOME=GENOME,
         ORGANISM=ORGANISM,
         ASSEMBLED_MOLECULES=ASSEMBLED_MOLECULES,
         CIRC_SEQS=CIRC_SEQS,
         GET_CHROM_SIZES=GET_CHROM_SIZES)
}

UCSC_registered_genomes <- function()
{
    dir_path <- system.file("registered_genomes", "UCSC",
                             package="GenomeInfoDb")
    file_paths <- list.files(dir_path, pattern="\\.R$", full.names=TRUE)
    genomes <- lapply(file_paths,
        function(file_path)
            as.list(.parse_script_for_UCSC_registered_genome(file_path)))
    colnames <- c("ORGANISM", "GENOME")
    listData <- lapply(setNames(colnames, tolower(colnames)),
        function(colname) {
            col <- vapply(genomes, `[[`, character(1), colname)
            if (colname == "ORGANISM")
                col <- factor(col)  # order of levels will dictate order
                                    # of rows in final DataFrame
            col
        }
    )
    listData$circ_seqs <- CharacterList(lapply(genomes, `[[`, "CIRC_SEQS"))
    ans <- S4Vectors:::new_DataFrame(listData, nrows=length(genomes))
    genome_trailing_digits <- sub("(.*[^0-9])([0-9]*)$", "\\2", ans$genome)
    oo <- order(ans$organism, as.integer(genome_trailing_digits))
    as.data.frame(ans[oo, , drop=FALSE])
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_chrom_info_from_UCSC()
###

.UCSC_cached_chrom_info <- new.env(parent=emptyenv())

.get_UCSC_chrom_info_for_unregistered_genome <- function(genome,
    assembled.molecules.only=FALSE,
    goldenPath.url=getOption("UCSC.goldenPath.url"),
    recache=FALSE)
{
    ans <- .UCSC_cached_chrom_info[[genome]]
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
        .UCSC_cached_chrom_info[[genome]] <- ans
    }
    if (assembled.molecules.only)
        warning(wmsg("'assembled.molecules' was ignored for unregistered ",
                     "genome ", genome, " (don't know what the assembled ",
                     "molecules are for unregistered genomes)"))
    ans
}

.get_UCSC_chrom_info_for_registered_genome <- function(genome, script_path,
    assembled.molecules.only=FALSE,
    goldenPath.url=getOption("UCSC.goldenPath.url"),
    recache=FALSE)
{
    vars <- .parse_script_for_UCSC_registered_genome(script_path)
    if (!identical(vars$GENOME, genome))
        stop(wmsg(script_path, ": script does not seem ",
                  "to be for genome ", genome))
    ASSEMBLED_MOLECULES <- vars$ASSEMBLED_MOLECULES

    ans <- .UCSC_cached_chrom_info[[genome]]
    if (is.null(ans) || recache) {
        GET_CHROM_SIZES <- vars$GET_CHROM_SIZES
        if (is.null(GET_CHROM_SIZES)) {
            ans <- fetch_chrom_sizes_from_UCSC(genome,
                                               goldenPath.url=goldenPath.url)
            stopifnot(nrow(ans) == length(ASSEMBLED_MOLECULES))
            oo <- match(ASSEMBLED_MOLECULES, ans[ , "chrom"])
            ans <- S4Vectors:::extract_data_frame_rows(ans, oo)
        } else {
            ans <- GET_CHROM_SIZES(goldenPath.url=goldenPath.url)
            stopifnot(is.data.frame(ans),
                      identical(sapply(ans, class),
                                c(chrom="character", size="integer")),
                      identical(ans[seq_along(ASSEMBLED_MOLECULES), "chrom"],
                                ASSEMBLED_MOLECULES))
        }
        ans$is_circular <- make_circ_flags_from_circ_seqs(ans[ , "chrom"],
                                                          vars$CIRC_SEQS)
        .UCSC_cached_chrom_info[[genome]] <- ans
    }
    if (assembled.molecules.only) {
        keep_idx <- seq_along(ASSEMBLED_MOLECULES)
        ans <- S4Vectors:::extract_data_frame_rows(ans, keep_idx)
    }
    ans
}

### Return a 3-column data.frame with columns "chrom" (character), "size"
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
    script_path <- system.file("registered_genomes", "UCSC", script_name,
                               package="GenomeInfoDb")
    if (identical(script_path, "")) {
        ans <- .get_UCSC_chrom_info_for_unregistered_genome(genome,
                    assembled.molecules.only=assembled.molecules.only,
                    goldenPath.url=goldenPath.url,
                    recache=recache)
    } else {
        ans <- .get_UCSC_chrom_info_for_registered_genome(genome, script_path,
                    assembled.molecules.only=assembled.molecules.only,
                    goldenPath.url=goldenPath.url,
                    recache=recache)
    }
    ans
}

