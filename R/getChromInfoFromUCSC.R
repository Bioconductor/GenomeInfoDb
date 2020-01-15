### =========================================================================
### getChromInfoFromUCSC()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### registered_UCSC_genomes()
###

.parse_script_for_registered_UCSC_genome <- function(script_path)
{
    filename <- basename(script_path)
    if (substr(filename, nchar(filename)-1L, nchar(filename)) != ".R")
        stop(wmsg("name of genome registration file '", filename, "' must ",
                  "have extension .R"))
    if (grepl(" ", filename, fixed=TRUE))
        stop(wmsg("name of genome registration file '", filename, "' must ",
                  "not contain spaces"))

    ## Placeholders. Will actually get defined when we source the script.
    ## See README.TXT in inst/registered_genomes/UCSC/ for the list of
    ## variables.
    GENOME <- ORGANISM <- ASSEMBLED_MOLECULES <- CIRC_SEQS <- NULL
    GET_CHROM_SIZES <- NCBI_LINKER <- NULL
    source(script_path, local=TRUE)

    stop_if <- function(notok, ...) {
        if (notok)
            stop("Error in UCSC genome registration file '", filename,
                 "':\n  ", wmsg(...))
    }

    ## Check GENOME.
    stop_if(is.null(GENOME), "'GENOME' must be defined")
    stop_if(!isSingleString(GENOME), "'GENOME' must be a single string")
    target <- substr(filename, 1L, nchar(filename)-2L)
    stop_if(!identical(target, GENOME), "'GENOME' must match filename")

    ## Check ORGANISM.
    stop_if(is.null(ORGANISM), "'ORGANISM' must be defined")
    stop_if(!isSingleString(ORGANISM), "'ORGANISM' must be a single string")
    stop_if(grepl("_", ORGANISM, fixed=TRUE),
            "underscores are not allowed in 'ORGANISM' (use spaces instead)")

    ## Check ASSEMBLED_MOLECULES.
    stop_if(is.null(ASSEMBLED_MOLECULES),
            "'ASSEMBLED_MOLECULES' must be defined")
    stop_if(!is.character(ASSEMBLED_MOLECULES),
            "'ASSEMBLED_MOLECULES' must be a character vector")
    stop_if(anyNA(ASSEMBLED_MOLECULES) ||
            !all(nzchar(ASSEMBLED_MOLECULES)) ||
            anyDuplicated(ASSEMBLED_MOLECULES),
            "'ASSEMBLED_MOLECULES' must ",
            "not contain NAs, empty strings, or duplicates")

    ## Check CIRC_SEQS.
    stop_if(is.null(CIRC_SEQS),
            "'CIRC_SEQS' must be defined")
    stop_if(!is.character(CIRC_SEQS),
            "'CIRC_SEQS' must be a character vector")
    stop_if(anyDuplicated(CIRC_SEQS),
            "'CIRC_SEQS' must not contain duplicates")
    stop_if(!all(CIRC_SEQS %in% ASSEMBLED_MOLECULES),
            "'CIRC_SEQS' must be a subset of 'ASSEMBLED_MOLECULES'")

    ## Check GET_CHROM_SIZES.
    if (!is.null(GET_CHROM_SIZES)) {
        stop_if(!is.function(GET_CHROM_SIZES),
                "when defined, 'GET_CHROM_SIZES' must be a function")
    }

    ## Check NCBI_LINKER.
    if (!is.null(NCBI_LINKER)) {
        stop_if(!is.list(NCBI_LINKER),
                "when defined, 'NCBI_LINKER' must be a named list")
        linker_fields <- names(NCBI_LINKER)
        stop_if(!is.character(linker_fields),
                "when defined, 'NCBI_LINKER' must be a named list")
        stop_if(anyNA(linker_fields) ||
                !all(nzchar(linker_fields)) ||
                anyDuplicated(linker_fields),
                "the names on 'NCBI_LINKER' must ",
                "not contain NAs, empty strings, or duplicates")
        stop_if(!("assembly_accession" %in% linker_fields),
                "'NCBI_LINKER' must have field \"assembly_accession\"")
        accession <- NCBI_LINKER$assembly_accession
        stop_if(!isSingleString(accession) || accession == "",
                "\"assembly_accession\" field in 'NCBI_LINKER' must ",
                "be a single non-empty string")
        assembly <- lookup_NCBI_accession2assembly(accession)
        stop_if(is.null(assembly),
                "\"assembly_accession\" field in 'NCBI_LINKER' must ",
                "be associated with a registered NCBI genome")
        stop_if(!identical(ORGANISM, assembly$organism),
                "the NCBI genome associated with the \"assembly_accession\" ",
                "field in 'NCBI_LINKER' is registered for an organism ",
                "(\"", assembly$organism, "\") that differs from 'ORGANISM' ",
                "(\"", ORGANISM, "\")")
    }

    list(GENOME=GENOME,
         ORGANISM=ORGANISM,
         ASSEMBLED_MOLECULES=ASSEMBLED_MOLECULES,
         CIRC_SEQS=CIRC_SEQS,
         GET_CHROM_SIZES=GET_CHROM_SIZES,
         NCBI_LINKER=NCBI_LINKER)
}

registered_UCSC_genomes <- function()
{
    dir_path <- system.file("registered_genomes", "UCSC",
                             package="GenomeInfoDb")
    file_paths <- list.files(dir_path, pattern="\\.R$", full.names=TRUE)
    assemblies <- lapply(file_paths, .parse_script_for_registered_UCSC_genome)

    colnames <- c("organism", "genome", "based_on_NCBI_genome", "circ_seqs")
    make_col <- function(j) {
        colname <- colnames[[j]]
        if (colname == "based_on_NCBI_genome") {
            col <- vapply(assemblies,
                function(assembly) {
                    linker <- assembly$NCBI_LINKER
                    if (is.null(linker))
                        return(NA_character_)
                    accession <- linker$assembly_accession
                    lookup_NCBI_accession2assembly(accession)$genome
                }, character(1), USE.NAMES=FALSE)
            return(col)
        }
        COLNAME <- toupper(colname)
        col0 <- lapply(assemblies, `[[`, COLNAME)
        if (colname == "circ_seqs")
            return(CharacterList(col0))
        stopifnot(all(lengths(col0) == 1L))
        col <- as.character(unlist(col0, use.names=FALSE))
        if (colname == "organism")
            col <- factor(col)  # order of levels will dictate order
                                # of rows in final DataFrame
        col
    }

    listData <- lapply(setNames(seq_along(colnames), colnames), make_col)
    ans <- S4Vectors:::new_DataFrame(listData, nrows=length(assemblies))
    genome_trailing_digits <- sub("(.*[^0-9])([0-9]*)$", "\\2", ans$genome)
    oo <- order(ans$organism, as.integer(genome_trailing_digits))
    as.data.frame(ans[oo, , drop=FALSE])
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getChromInfoFromUCSC()
###

.UCSC_cached_chrom_info <- new.env(parent=emptyenv())

.get_chrom_info_for_unregistered_UCSC_genome <- function(genome,
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
            stop(wmsg("unknown UCSC genome assembly: ", genome))
        oo <- orderSeqlevels(ans[ , "chrom"])
        ans <- S4Vectors:::extract_data_frame_rows(ans, oo)

        ## Add columns "assembled" and "circular".
        assembled <- rep.int(NA, nrow(ans))
        circular <- make_circ_flags_from_circ_seqs(ans[ , "chrom"])
        ans <- cbind(ans, assembled=assembled, circular=circular)

        .UCSC_cached_chrom_info[[genome]] <- ans
    }
    if (assembled.molecules.only)
        warning(wmsg("'assembled.molecules' was ignored for unregistered ",
                     "genome ", genome, " (don't know what the assembled ",
                     "molecules are for unregistered UCSC genomes)"))
    ans
}

.get_chrom_info_for_registered_UCSC_genome <- function(script_path,
    assembled.molecules.only=FALSE,
    goldenPath.url=getOption("UCSC.goldenPath.url"),
    recache=FALSE)
{
    vars <- .parse_script_for_registered_UCSC_genome(script_path)
    GENOME <- vars$GENOME
    ASSEMBLED_MOLECULES <- vars$ASSEMBLED_MOLECULES
    nb_assembled <- length(ASSEMBLED_MOLECULES)
    assembled_idx <- seq_len(nb_assembled)

    ans <- .UCSC_cached_chrom_info[[GENOME]]
    if (is.null(ans) || recache) {
        GET_CHROM_SIZES <- vars$GET_CHROM_SIZES
        if (is.null(GET_CHROM_SIZES)) {
            ans <- fetch_chrom_sizes_from_UCSC(GENOME,
                                               goldenPath.url=goldenPath.url)
            stopifnot(nrow(ans) == nb_assembled)
            oo <- match(ASSEMBLED_MOLECULES, ans[ , "chrom"])
            stopifnot(!anyNA(oo))
            ans <- S4Vectors:::extract_data_frame_rows(ans, oo)
        } else {
            ans <- GET_CHROM_SIZES(goldenPath.url=goldenPath.url)
            stopifnot(is.data.frame(ans),
                      identical(sapply(ans, class),
                                c(chrom="character", size="integer")),
                      identical(head(ans[ , "chrom"], n=nb_assembled),
                                ASSEMBLED_MOLECULES),
                      !anyNA(ans[ , "chrom"]),
                      all(nzchar(ans[ , "chrom"])),
                      !anyDuplicated(ans[ , "chrom"]))
        }

        ## Add columns "assembled" and "circular".
        assembled <- logical(nrow(ans))
        assembled[assembled_idx] <- TRUE
        circular <- make_circ_flags_from_circ_seqs(ans[ , "chrom"],
                                                   vars$CIRC_SEQS)
        ans <- cbind(ans, assembled=assembled, circular=circular)

        .UCSC_cached_chrom_info[[GENOME]] <- ans
    }
    if (assembled.molecules.only)
        ans <- S4Vectors:::extract_data_frame_rows(ans, assembled_idx)
    ans
}

### Return a 4-column data.frame with columns "chrom" (character), "size"
### (integer), "assembled" (logical), and "circular" (logical).
getChromInfoFromUCSC <- function(genome,
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
        ans <- .get_chrom_info_for_unregistered_UCSC_genome(genome,
                    assembled.molecules.only=assembled.molecules.only,
                    goldenPath.url=goldenPath.url,
                    recache=recache)
    } else {
        ans <- .get_chrom_info_for_registered_UCSC_genome(script_path,
                    assembled.molecules.only=assembled.molecules.only,
                    goldenPath.url=goldenPath.url,
                    recache=recache)
    }
    ans
}

