### =========================================================================
### getChromInfoFromUCSC()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .add_NCBI_columns()
###

.safe_match <- function(query, NCBI_accns)
{
    hits <- findMatches(query, NCBI_accns)
    q_hits <- queryHits(hits)
    s_hits <- subjectHits(hits)
    ambig_q_idx <- which(duplicated(q_hits))
    if (length(ambig_q_idx) != 0L) {
        ambig_idx <- unique(q_hits[ambig_q_idx])
        in1string <- paste0(names(query)[ambig_idx], collapse=", ")
        stop(wmsg("UCSC seqlevel(s) matching more than one accession ",
                  "number: ", in1string))
    }
    ambig_s_idx <- which(duplicated(s_hits))
    if (length(ambig_s_idx) != 0L) {
        ambig_idx <- unique(s_hits[ambig_s_idx])
        in1string <- paste0(NCBI_accns[ambig_idx], collapse=", ")
        stop(wmsg("accession number(s) with more than one UCSC seqlevel ",
                  "match: ", in1string))
    }
    ans <- rep.int(NA_integer_, length(query))
    ans[q_hits] <- s_hits
    ans
}

.match_UCSC_seqlevels_to_NCBI_accns <- function(UCSC_seqlevels, NCBI_accns,
                                                accn_suffix="")
{
    query <- chartr("-", ".", UCSC_seqlevels)
    query <- paste0(query, accn_suffix)
    names(query) <- UCSC_seqlevels
    .safe_match(tolower(query), tolower(NCBI_accns))
}

.match_trimmed_UCSC_random_seqlevels_to_NCBI_accns <-
    function(UCSC_seqlevels, NCBI_accns)
{
    query <- sub("^[^_]*_", "", sub("_random$", "", UCSC_seqlevels))
    query <- chartr("v", ".", query)
    names(query) <- UCSC_seqlevels
    .safe_match(tolower(query), tolower(NCBI_accns))
}

.match_UCSC_seqlevels_part2_to_NCBI_accns <-
    function(UCSC_seqlevels, NCBI_accns, accn_prefix="")
{
    ans <- rep.int(NA_integer_, length(UCSC_seqlevels))
    seqlevel_parts <- strsplit(UCSC_seqlevels, "_")
    nparts <- lengths(seqlevel_parts)
    idx2 <- which(nparts >= 2L)
    if (length(idx2) == 0L)
        return(ans)
    offsets <- c(0L, cumsum(nparts[idx2[-length(idx2)]]))
    query <- unlist(seqlevel_parts[idx2], use.names=FALSE)[offsets + 2L]
    query <- chartr("v", ".", query)
    unversioned_idx <- grep(".", query, fixed=TRUE, invert=TRUE)
    if (length(unversioned_idx) != 0L)
        query[unversioned_idx] <- paste0(query[unversioned_idx], ".1")
    query <- paste0(accn_prefix, query)
    names(query) <- UCSC_seqlevels[idx2]
    ans[idx2] <- .safe_match(tolower(query), tolower(NCBI_accns))
    ans
}

### 'NCBI_seqlevels', 'NCBI_gbkaccns', and 'NCBI_refseqaccns' must be
### parallel vectors.
.map_UCSC_seqlevels_to_NCBI_seqlevels <- function(UCSC_seqlevels,
                                                  NCBI_seqlevels,
                                                  NCBI_ucsc_names,
                                                  NCBI_gbkaccns,
                                                  NCBI_refseqaccns,
                                                  special_mappings=NULL)
{
    ans <- rep.int(NA_integer_, length(UCSC_seqlevels))

    ## 1. Handle special mappings.
    if (!is.null(special_mappings)) {
        m1 <- match(names(special_mappings), UCSC_seqlevels)
        if (any(is.na(m1)))
            stop(wmsg("'special_mappings' contains sequence names ",
                      "not in 'UCSC_seqlevels'"))
        m2 <- match(special_mappings, NCBI_seqlevels)
        if (any(is.na(m2)))
            stop(wmsg("'special_mappings' has values not in 'NCBI_seqlevels'"))
        ans[m1] <- m2
    }
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 2. We assign based on NCBI_ucsc_names.
    m <- match(UCSC_seqlevels[unmapped_idx], NCBI_ucsc_names)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 3. We assign based on exact matching (case insensitive) of the
    ##    seqlevels.
    ucsc_seqlevels <- tolower(UCSC_seqlevels)
    ncbi_seqlevels <- tolower(NCBI_seqlevels)
    m <- match(ucsc_seqlevels[unmapped_idx], ncbi_seqlevels)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 4. We assign based on exact matching (case insensitive) of the
    ##    seqlevels after removal of the "chr" prefix.
    nochr_ucsc_seqlevels <- sub("^chr", "", ucsc_seqlevels[unmapped_idx])
    nochr_ncbi_seqlevels <- sub("^chr", "", ncbi_seqlevels)
    m <- match(nochr_ucsc_seqlevels, nochr_ncbi_seqlevels)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 5. We assign based on GenBank accession number found in UCSC seqlevels.
    m <- .match_UCSC_seqlevels_to_NCBI_accns(UCSC_seqlevels[unmapped_idx],
                                             NCBI_gbkaccns)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 6. We assign based on GenBank accession number found in UCSC seqlevels
    ##    after adding .1 suffix to it.
    m <- .match_UCSC_seqlevels_to_NCBI_accns(UCSC_seqlevels[unmapped_idx],
                                             NCBI_gbkaccns, accn_suffix=".1")
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 7. We assign based on RefSeq accession number found in UCSC seqlevels
    ##    after trimming the first part (e.g. "chr1_") and "_random" suffix
    ##    from the seqlevels.
    m <- .match_trimmed_UCSC_random_seqlevels_to_NCBI_accns(
                                                 UCSC_seqlevels[unmapped_idx],
                                                 NCBI_refseqaccns)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 8. We assign based on GenBank accession number found in part 2 of
    ##    UCSC seqlevels.
    m <- .match_UCSC_seqlevels_part2_to_NCBI_accns(
                                        UCSC_seqlevels[unmapped_idx],
                                        NCBI_gbkaccns)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 9. We assign based on GenBank accession number found in part 2 of
    ##    UCSC seqlevels after adding AAD prefix to it.
    ans[unmapped_idx] <- .match_UCSC_seqlevels_part2_to_NCBI_accns(
                                    UCSC_seqlevels[unmapped_idx],
                                    NCBI_gbkaccns,
                                    accn_prefix="AAD")
    ans
}

.add_NCBI_columns <- function(UCSC_chrom_info, assembly_accession,
    AssemblyUnits=NULL, special_mappings=NULL, unmapped_seqs=NULL,
    drop_unmapped=FALSE)
{
    UCSC_seqlevels <- UCSC_chrom_info[ , "chrom"]

    if (length(unmapped_seqs) != 0L) {
        unmapped_seqs_role <- rep.int(names(unmapped_seqs),
                                      lengths(unmapped_seqs))
        unmapped_seqs <- unlist(unmapped_seqs, use.names=FALSE)
        unmapped_idx <- match(unmapped_seqs, UCSC_seqlevels)
        stopifnot(!any(is.na(unmapped_idx)))
    }

    NCBI_chrom_info <- getChromInfoFromNCBI(assembly_accession,
                                            assembly.units=AssemblyUnits)

    NCBI_ucsc_names <- NCBI_chrom_info[ , "UCSCStyleName"]
    NCBI_seqlevels <- NCBI_chrom_info[ , "SequenceName"]
    NCBI_gbkaccns <- NCBI_chrom_info[ , "GenBankAccn"]
    NCBI_refseqaccns <- NCBI_chrom_info[ , "RefSeqAccn"]
    oo <- .map_UCSC_seqlevels_to_NCBI_seqlevels(UCSC_seqlevels, NCBI_seqlevels,
                                      NCBI_ucsc_names,
                                      NCBI_gbkaccns, NCBI_refseqaccns,
                                      special_mappings=special_mappings)
    oo_is_NA <- is.na(oo)
    mapped_idx <- which(!oo_is_NA)

    if (length(unmapped_seqs) != 0L)
        stopifnot(all(oo_is_NA[unmapped_idx]))

    if (isTRUE(drop_unmapped)) {
        UCSC_chrom_info <-
            S4Vectors:::extract_data_frame_rows(UCSC_chrom_info, mapped_idx)
        oo <- oo[mapped_idx]
        mapped_idx <- seq_along(oo)
        ## Before we drop the unmapped UCSC seqlevels, we want to make sure
        ## that all the NCBI seqlevels are reverse mapped. For example, in
        ## the case of hg38, the chromInfo table at UCSC contains sequences
        ## that belong to GRCh38.p12 but not to GRCh38. So we want to make
        ## sure that after we drop these "foreign sequences", we are left
        ## with a one-to-one mapping between UCSC seqlevels and the 455 NCBI
        ## seqlevels that are in the assembly report for GRCh38.
        idx <- setdiff(seq_along(NCBI_seqlevels), oo)
        if (length(idx) != 0L) {
            in1string <- paste0(NCBI_seqlevels[idx], collapse=", ")
            stop(wmsg("no UCSC seqlevel could be mapped to the following ",
                      "NCBI seqlevel(s): ", in1string))
        }
    } else {
        unexpectedly_unmapped_idx <-
            which(oo_is_NA & !(UCSC_seqlevels %in% unmapped_seqs))
        if (length(unexpectedly_unmapped_idx) != 0L) {
            in1string <- paste0(UCSC_seqlevels[unexpectedly_unmapped_idx],
                                collapse=", ")
            stop(wmsg("cannot map the following UCSC seqlevel(s) to an ",
                      "NCBI seqlevel: ", in1string))
        }
    }

    NCBI_chrom_info <- S4Vectors:::extract_data_frame_rows(NCBI_chrom_info, oo)

    ## NCBI columns "circular", "SequenceLength", and "UCSCStyleName"
    ## are expected to be redundant with UCSC columns "circular", "size",
    ## and "chrom", respectively. Let's make sure these 3 NCBI columns
    ## agree with their UCSC counterpart before dropping them.

    stopifnot(identical(UCSC_chrom_info[mapped_idx, "circular"],
                        NCBI_chrom_info[mapped_idx, "circular"]))

    compare_idx <- which(!is.na(NCBI_chrom_info[ , "SequenceLength"]))
    stopifnot(identical(UCSC_chrom_info[compare_idx, "size"],
                        NCBI_chrom_info[compare_idx, "SequenceLength"]))

    compare_idx <- which(!is.na(NCBI_chrom_info[ , "UCSCStyleName"]))
    stopifnot(identical(UCSC_chrom_info[compare_idx, "chrom"],
                        NCBI_chrom_info[compare_idx, "UCSCStyleName"]))

    keep_NCBI_colnames <- c("SequenceName", "SequenceRole",
                            "GenBankAccn", "Relationship", "RefSeqAccn",
                            "AssemblyUnit")
    ans <- cbind(UCSC_chrom_info, NCBI_chrom_info[ , keep_NCBI_colnames])

    if (length(unmapped_seqs) != 0L)
        ans[unmapped_idx, "SequenceRole"] <- unmapped_seqs_role

    stopifnot(!is.unsorted(ans[ , "SequenceRole"]))
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .add_ensembl_column()
###

.fetch_ucscToEnsembl_table <- function(genome,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    col2class <- c(ucsc="character", ensembl="character")
    fetch_table_from_UCSC(genome, "ucscToEnsembl",
                          col2class=col2class,
                          goldenPath.url=goldenPath.url)
}

### Filters and reformats the chromAlias data to so it's returned in the
### same format as the ucscToEnsembl table (i.e. 2-column data frame with
### columns "ucsc" and "ensembl").
.fetch_ucsc2ensembl_from_chromAlias_table <- function(genome,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    col2class <- c(ensembl="character", ucsc="character", source="factor")
    chromAlias <- fetch_table_from_UCSC(genome, "chromAlias",
                                        col2class=col2class,
                                        goldenPath.url=goldenPath.url)
    ## Filters and reformats.
    keep_idx <- grep("ensembl", chromAlias[ , "source"], fixed=TRUE)
    ucsc2ensembl <- chromAlias[keep_idx, c("ucsc", "ensembl")]

    ## In some **very rare** situations, it seems that the chromAlias table
    ## can contain ambiguous ucsc-to-ensembl mappings where a given UCSC
    ## chromosome name is mapped to more than one Ensembl name. For example
    ## this happens for rn6 where unplaced scaffold chrUn_AABR07022993v1 is
    ## mapped to Ensembl names AABR07022518.1, AABR07023006.1, AABR07022993.1.
    ## Note that as of Jan 2020, the chromAlias table for rn6 is the only
    ## known case of ambiguous ucsc-to-ensembl mapping.
    if (anyDuplicated(ucsc2ensembl[ , "ucsc"])) {
        ## The disambiguation algorithm below assumes that the ambiguous
        ## ucsc-to-ensembl mappings are of the form xxYYY-to-ZZZ where
        ## ZZZ is a GenBank accession, YYY a GenBank accession where the
        ## dot (".") has been replaced with a "v" and xx can be any suffix.
        ## To disambiguate, we only keep mappings where YYY and ZZZ are the
        ## same (modulo the dot/v substitution).
        idx <- which(duplicated(ucsc2ensembl[ , "ucsc"]) |
                     duplicated(ucsc2ensembl[ , "ucsc"], fromLast=TRUE))
        ensembl <- chartr(".", "v", ucsc2ensembl[idx , "ensembl"])
        ucsc <- ucsc2ensembl[idx , "ucsc"]
        ucsc_nc <- nchar(ucsc)
        ucsc_suffix <- substr(ucsc, ucsc_nc - nchar(ensembl) + 1L, ucsc_nc)
        ok <- ucsc_suffix == ensembl
        if (sum(ok) == 0L)
            stop(wmsg("failed to disambiguate 'ucsc2ensembl'"))
        if (!all(ucsc %in% ucsc[ok]))
            stop(wmsg("failed to disambiguate 'ucsc2ensembl'"))
        if (anyDuplicated(ucsc[ok]))
            stop(wmsg("failed to disambiguate 'ucsc2ensembl'"))
        drop_idx <- idx[!ok]
        ucsc2ensembl <- ucsc2ensembl[-drop_idx, , drop=FALSE]
    }
    ucsc2ensembl
}

.try_to_fetch_ucscToEnsembl <- function(genome,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    message("Trying 'ucscToEnsembl' ... ", appendLF=FALSE)
    ucscToEnsembl <- try(
        suppressWarnings(
            .fetch_ucscToEnsembl_table(genome, goldenPath.url=goldenPath.url)
        ), silent=TRUE)
    if (inherits(ucscToEnsembl, "try-error")) {
        message("FAILED! (table does not exist)")
        return(NULL)
    }
    message("OK")
    ucscToEnsembl
}

.try_to_fetch_ucsc2ensembl_from_chromAlias <- function(genome,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    message("Trying 'chromAlias' ... ", appendLF=FALSE)
    ucsc2ensembl <- try(
        suppressWarnings(
            .fetch_ucsc2ensembl_from_chromAlias_table(genome,
                                                goldenPath.url=goldenPath.url)
        ), silent=TRUE)
    if (inherits(ucsc2ensembl, "try-error")) {
        message("FAILED! (table does not exist)")
        return(NULL)
    }
    if (nrow(ucsc2ensembl) == 0L) {
        message("FAILED! (table exists but has no Ensembl information)")
        return(NULL)
    }
    message("OK")
    ucsc2ensembl
}

.fetch_Ensembl_column <- function(UCSC_chroms, genome, table=NA,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    if (identical(table, "ucscToEnsembl")) {
        ucsc2ensembl <- .fetch_ucscToEnsembl_table(genome,
                                                goldenPath.url=goldenPath.url)
    } else if (identical(table, "chromAlias")) {
        ucsc2ensembl <- .fetch_ucsc2ensembl_from_chromAlias_table(genome,
                                                goldenPath.url=goldenPath.url)
    } else if (identical(table, NA)) {
        table <- "ucscToEnsembl"
        ucsc2ensembl <- .try_to_fetch_ucscToEnsembl(genome,
                                      goldenPath.url=goldenPath.url)
        if (is.null(ucsc2ensembl)) {
            table <- "chromAlias"
            ucsc2ensembl <- .try_to_fetch_ucsc2ensembl_from_chromAlias(genome,
                                      goldenPath.url=goldenPath.url)
            if (is.null(ucsc2ensembl)) {
                warning(wmsg("No Ensembl sequence names could be found ",
                             "for UCSC genome ", genome))
                return(rep.int(NA_character_, length(UCSC_chroms)))
            }
        }
    } else {
        stop(wmsg("supplied 'table' is not supported"))
    }
    oo <- match(UCSC_chroms, ucsc2ensembl[ , "ucsc"])
    if (is.null(table) && anyNA(oo))
        warning(wmsg("incomplete data in '", table, "' table for ",
                     "UCSC genome ", genome, ": some UCSC sequence ",
                     "names are missing"))
    ucsc2ensembl[oo , "ensembl"]
}

.add_ensembl_column <- function(UCSC_chrom_info, genome, table=NA,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    ensembl <- .fetch_Ensembl_column(UCSC_chrom_info[ , "chrom"],
                                     genome, table=table,
                                     goldenPath.url=goldenPath.url)
    cbind(UCSC_chrom_info, ensembl=ensembl, stringsAsFactors=FALSE)
}


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
    GET_CHROM_SIZES <- NCBI_LINKER <- ENSEMBL_LINKER <- NULL
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
        NCBI_assembly <- lookup_NCBI_accession2assembly(accession)
        stop_if(is.null(NCBI_assembly),
                "\"assembly_accession\" field in 'NCBI_LINKER' must ",
                "be associated with a registered NCBI genome")
        stop_if(!identical(ORGANISM, NCBI_assembly$organism),
                "the NCBI genome associated with the \"assembly_accession\" ",
                "field in 'NCBI_LINKER' is registered for an organism ",
                "(\"", NCBI_assembly$organism, "\") that differs ",
                "from 'ORGANISM' (\"", ORGANISM, "\")")
    }

    ## Check ENSEMBL_LINKER.
    if (!is.null(ENSEMBL_LINKER)) {
        stop_if(!isSingleString(ENSEMBL_LINKER) || ENSEMBL_LINKER == "",
                "when defined, 'ENSEMBL_LINKER' must ",
                "be a single non-empty string")
    }

    list(GENOME=GENOME,
         ORGANISM=ORGANISM,
         ASSEMBLED_MOLECULES=ASSEMBLED_MOLECULES,
         CIRC_SEQS=CIRC_SEQS,
         GET_CHROM_SIZES=GET_CHROM_SIZES,
         NCBI_LINKER=NCBI_LINKER,
         ENSEMBL_LINKER=ENSEMBL_LINKER)
}

registered_UCSC_genomes <- function()
{
    dir_path <- system.file("registered_genomes", "UCSC",
                             package="GenomeInfoDb")
    file_paths <- list.files(dir_path, pattern="\\.R$", full.names=TRUE)
    assemblies <- lapply(file_paths, .parse_script_for_registered_UCSC_genome)

    colnames <- c("organism", "genome", "based_on_NCBI_genome",
                  "assembly_accession", "with_Ensembl", "circ_seqs")
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
        if (colname == "assembly_accession") {
            col <- vapply(assemblies,
                function(assembly) {
                    linker <- assembly$NCBI_LINKER
                    if (is.null(linker))
                        return(NA_character_)
                    linker$assembly_accession
                }, character(1), USE.NAMES=FALSE)
            return(col)
        }
        if (colname == "with_Ensembl") {
            col <- vapply(assemblies,
                function(assembly) {
                    linker <- assembly$ENSEMBL_LINKER
                    !is.null(linker)
                }, logical(1), USE.NAMES=FALSE)
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
    add.ensembl.col=FALSE,
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
    if (isTRUE(add.ensembl.col)) {
        warning(wmsg("'add.ensembl.col' got ignored for unregistered ",
                     "UCSC genome ", genome, "."),
                     "\n  ",
                wmsg("Use 'add.ensembl.col=\"force\"' to try adding ",
                     "the \"ensembl\" column anyway."))
    } else if (identical(add.ensembl.col, "force")) {
        ans <- .add_ensembl_column(ans, genome,
                                   goldenPath.url=goldenPath.url)
    }
    if (assembled.molecules.only)
        warning(wmsg("'assembled.molecules' got ignored for unregistered ",
                     "UCSC genome ", genome, " (don't know what the ",
                     "assembled molecules are for an unregistered UCSC ",
                     "genome)"))
    ans
}

.get_NCBI_linker <- function(add.NCBI.cols, linker, GENOME)
{
    if (is.list(add.NCBI.cols))
        return(add.NCBI.cols)
    if (isFALSE(add.NCBI.cols))
        return(NULL)
    if (is.null(linker))
        warning(wmsg("'add.NCBI.cols' got ignored for registered ",
                     "UCSC genome ", GENOME, " (additional NCBI columns ",
                     "are only available for registered UCSC genomes ",
                     "based on an NCBI genome)"))
    linker
}

.get_Ensembl_linker <- function(add.ensembl.col, linker, GENOME)
{
    if (identical(add.ensembl.col, "force"))
        return(NA)
    if (isFALSE(add.ensembl.col))
        return(NULL)
    if (is.null(linker))
        warning(wmsg("'add.ensembl.col' got ignored for registered ",
                     "UCSC genome ", GENOME, " (genome was registered ",
                     "with no Ensembl sequence names linked to it)."),
                     "\n  ",
                wmsg("Use 'add.ensembl.col=\"force\"' to try adding ",
                     "the \"ensembl\" column anyway."))
    linker
}

.get_chrom_info_for_registered_UCSC_genome <- function(script_path,
    assembled.molecules.only=FALSE,
    add.NCBI.cols=FALSE,
    add.ensembl.col=FALSE,
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

    ## Add NCBI cols.
    NCBI_linker <- .get_NCBI_linker(add.NCBI.cols, vars$NCBI_LINKER, GENOME)
    if (!is.null(NCBI_linker))
        ans <- do.call(.add_NCBI_columns, c(list(ans), NCBI_linker))

    ## Add Ensembl col.
    Ensemb_linker <- .get_Ensembl_linker(add.ensembl.col, vars$ENSEMBL_LINKER,
                                         GENOME)
    if (!is.null(Ensemb_linker))
        ans <- .add_ensembl_column(ans, GENOME, table=Ensemb_linker,
                                   goldenPath.url=goldenPath.url)

    if (assembled.molecules.only)
        ans <- S4Vectors:::extract_data_frame_rows(ans, assembled_idx)
    ans
}

### Return a 4-column data.frame with columns "chrom" (character), "size"
### (integer), "assembled" (logical), and "circular" (logical).
getChromInfoFromUCSC <- function(genome,
    assembled.molecules.only=FALSE,
    add.NCBI.cols=FALSE,
    add.ensembl.col=FALSE,
    goldenPath.url=getOption("UCSC.goldenPath.url"),
    recache=FALSE,
    as.Seqinfo=FALSE)
{
    if (!isSingleString(genome))
        stop(wmsg("'genome' must be a single string"))
    if (!isTRUEorFALSE(assembled.molecules.only))
        stop(wmsg("'assembled.molecules.only' must be TRUE or FALSE"))
    if (!(isTRUEorFALSE(add.NCBI.cols) ||
          is.list(add.NCBI.cols) && !is.null(names(add.NCBI.cols))))
        stop(wmsg("'add.NCBI.cols' must be TRUE or FALSE ",
                  "(experts can also supply an \"NCBI linker\" ",
                  "as a named list)"))
    if (!(isTRUEorFALSE(add.ensembl.col) ||
          identical(add.ensembl.col, "force")))
        stop(wmsg("'add.ensembl.col' must be TRUE or FALSE ",
                  "(experts can also use \"force\")"))
    if (!isSingleString(goldenPath.url))
        stop(wmsg("'goldenPath.url' must be a single string"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(as.Seqinfo))
        stop(wmsg("'as.Seqinfo' must be TRUE or FALSE"))

    script_name <- paste0(genome, ".R")
    script_path <- system.file("registered_genomes", "UCSC", script_name,
                               package="GenomeInfoDb")
    if (identical(script_path, "")) {
        if (!isFALSE(add.NCBI.cols))
            warning(wmsg("'add.NCBI.cols' got ignored for unregistered ",
                         "UCSC genome ", genome, " (additional NCBI columns ",
                         "are only available for registered UCSC genomes ",
                         "based on an NCBI genome)"))
        ans <- .get_chrom_info_for_unregistered_UCSC_genome(genome,
                    assembled.molecules.only=assembled.molecules.only,
                    add.ensembl.col=add.ensembl.col,
                    goldenPath.url=goldenPath.url,
                    recache=recache)
    } else {
        ans <- .get_chrom_info_for_registered_UCSC_genome(script_path,
                    assembled.molecules.only=assembled.molecules.only,
                    add.NCBI.cols=add.NCBI.cols,
                    add.ensembl.col=add.ensembl.col,
                    goldenPath.url=goldenPath.url,
                    recache=recache)
    }
    if (!as.Seqinfo)
        return(ans)
    Seqinfo(seqnames=ans[ , "chrom"],
            seqlengths=ans[ , "size"],
            isCircular=ans[ , "circular"],
            genome=genome)
}

