### =========================================================================
### getChromInfoFromNCBI()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load and access db of registered NCBI genomes
###

.NCBI_genome2accession <- new.env(parent=emptyenv())
.NCBI_accession2assembly <- new.env(parent=emptyenv())

.load_registered_NCBI_genome <- function(file_path)
{
    filename <- basename(file_path)
    if (substr(filename, nchar(filename)-1L, nchar(filename)) != ".R")
        stop(wmsg("name of genome registration file '", filename, "' must ",
                  "have extension .R"))
    if (grepl(" ", filename, fixed=TRUE))
        stop(wmsg("name of genome registration file '", filename, "' must ",
                  "not contain spaces"))

    ## Placeholders. Will actually get defined when we source the
    ## assembly files.
    ORGANISM <- NULL    # Expected to be a single string.
    ASSEMBLIES <- NULL  # Expected to be a list with one list element per
                        # assembly.
    source(file_path, local=TRUE)

    stop_if <- function(notok, ...) {
        if (notok)
            stop("Error in NCBI genome registration file '", filename,
                 "':\n  ", wmsg(...))
    }

    ## Check ORGANISM.
    stop_if(is.null(ORGANISM), "'ORGANISM' must be defined")
    stop_if(!isSingleString(ORGANISM), "'ORGANISM' must be a single string")
    stop_if(grepl("_", ORGANISM, fixed=TRUE),
            "underscores are not allowed in 'ORGANISM' (use spaces instead)")
    target <- chartr("_", " ", substr(filename, 1L, nchar(filename)-2L))
    stop_if(!identical(target, ORGANISM),
            "'ORGANISM' must match filename ",
            "(with spaces in place of underscores)")

    ## Check ASSEMBLIES.
    stop_if(is.null(ASSEMBLIES), "'ASSEMBLIES' must be defined")
    stop_if(!is.list(ASSEMBLIES), "'ASSEMBLIES' must be a list of lists")

    required_fields <- c("genome", "date", "assembly_accession", "circ_seqs")
    for (i in seq_along(ASSEMBLIES)) {
        assembly <- ASSEMBLIES[[i]]

        stop_if(!is.list(assembly),
                "'ASSEMBLIES[[", i, "]]' must be a named list")
        assembly_fields <- names(assembly)
        stop_if(!is.character(assembly_fields),
                "'ASSEMBLIES[[", i, "]]' must be a named list")
        stop_if(anyNA(assembly_fields) ||
                !all(nzchar(assembly_fields)) ||
                anyDuplicated(assembly_fields),
                "the names on 'ASSEMBLIES[[", i, "]]' must ",
                "not contain NAs, empty strings, or duplicates")
        stop_if(!all(required_fields %in% names(assembly)),
                "'ASSEMBLIES[[", i, "]]' must have fields: ",
                paste(paste0("\"", required_fields, "\""), collapse=", "))

        ## Check "genome" field (required).
        genome <- assembly$genome
        stop_if(!isSingleString(genome) || genome == "",
                "\"genome\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "be a single non-empty string")

        ## Check "date" field (required).
        date <- assembly$date
        stop_if(!isSingleString(date) || date == "",
                "\"date\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "be a single non-empty string")

        ## Check "accession" field (required).
        accession <- assembly$assembly_accession
        stop_if(!isSingleString(accession) || accession == "",
                "\"accession\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "be a single non-empty string")
        stop_if(!is.null(.NCBI_accession2assembly[[accession]]),
                "assembly accession ", accession, " used more than once")

        ## Check "circ_seqs" field (required).
        circ_seqs <- assembly$circ_seqs
        stop_if(!is.character(circ_seqs),
                "\"circ_seqs\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "be a character vector")
        stop_if(anyNA(circ_seqs) ||
                !all(nzchar(circ_seqs)) ||
                anyDuplicated(circ_seqs),
                "\"circ_seqs\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "not contain NAs, empty strings, or duplicates")

        ## Check optional fields.
        extra_info <- assembly$extra_info
        if (!is.null(extra_info)) {
            stop_if(!is.character(extra_info),
                "\"extra_info\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "be a named character vector")
            stop_if(anyNA(extra_info) || !all(nzchar(extra_info)),
                "\"extra_info\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "not contain NAs or empty strings")
            tags <- names(extra_info)
            stop_if(!is.character(tags),
                "\"extra_info\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "be a named character vector")
            stop_if(anyNA(tags) || anyDuplicated(tags),
                "the names on \"extra_info\" field 'ASSEMBLIES[[", i, "]]' ",
                "must not contain NAs or duplicates")
        }

        genome <- tolower(genome)
        .NCBI_genome2accession[[genome]] <-
            c(.NCBI_genome2accession[[genome]], accession)
        assembly$organism <- ORGANISM
        assembly$rank_within_organism <- i
        .NCBI_accession2assembly[[accession]] <- assembly
    }
}

.load_registered_NCBI_genomes <- function()
{
    dir_path <- system.file("registered_genomes", "NCBI",
                             package="GenomeInfoDb")
    file_paths <- list.files(dir_path, pattern="\\.R$", full.names=TRUE)
    for (file_path in file_paths)
        .load_registered_NCBI_genome(file_path)
}

.make_extra_info_col <- function(col0, ifmissing="", sep=":", collapse="|")
{
    stopifnot(is.list(col0), is.character(ifmissing), length(ifmissing) == 1L)
    col <- rep.int(ifmissing, length(col0))
    idx1 <- which(lengths(col0) != 0L)
    if (length(idx1) != 0L) {
        col1 <- col0[idx1]
        unlisted_col1 <- unlist(col1)
        stopifnot(is.character(unlisted_col1))
        tags <- names(unlisted_col1)
        stopifnot(is.character(tags),
                  !anyNA(tags),
                  all(nzchar(tags)))
        col1 <- relist(paste(tags, unlisted_col1, sep=sep), col1)
        col[idx1] <- unstrsplit(col1, sep=collapse)
    }
    col
}

registered_NCBI_genomes <- function()
{
    if (length(.NCBI_accession2assembly) == 0L)
        .load_registered_NCBI_genomes()
    assemblies <- unname(as.list(.NCBI_accession2assembly, all.names=TRUE))

    colnames <- c("organism", "rank_within_organism", "genome", "date",
                  "extra_info", "assembly_accession", "circ_seqs")
    make_col <- function(colname) {
        col0 <- lapply(assemblies, `[[`, colname)
        if (colname == "extra_info")
            return(factor(.make_extra_info_col(col0)))
        if (colname == "circ_seqs")
            return(CharacterList(col0))
        stopifnot(all(lengths(col0) == 1L))
        col <- unlist(col0, use.names=FALSE)
        if (colname == "rank_within_organism")
            return(as.integer(col))
        col <- as.character(col)
        if (colname == "organism")
            col <- factor(col)  # order of levels will dictate order
                                # of rows in final DataFrame
        col
    }

    listData <- lapply(setNames(colnames, colnames), make_col)
    ans <- S4Vectors:::new_DataFrame(listData, nrows=length(assemblies))
    oo <- order(ans$organism, ans$rank_within_organism)
    as.data.frame(drop_cols(ans, "rank_within_organism")[oo, , drop=FALSE])
}

lookup_NCBI_genome2accession <- function(genome)
{
    stopifnot(isSingleString(genome))
    if (length(.NCBI_genome2accession) == 0L)
        .load_registered_NCBI_genomes()
    .NCBI_genome2accession[[tolower(genome)]]
}

lookup_NCBI_accession2assembly <- function(accession)
{
    stopifnot(isSingleString(accession))
    if (length(.NCBI_accession2assembly) == 0L)
        .load_registered_NCBI_genomes()
    .NCBI_accession2assembly[[accession]]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getChromInfoFromNCBI()
###

.format_assembly_report <- function(assembly_report, circ_seqs=NULL)
{
    drop_columns <- c("AssignedMolecule", "AssignedMoleculeLocationOrType")
    ans <- assembly_report[ , !(colnames(assembly_report) %in% drop_columns)]
    ans[ , "SequenceName"] <- as.character(ans[ , "SequenceName"])
    SequenceRole_levels <- c("assembled-molecule",
                             "alt-scaffold",
                             "unlocalized-scaffold",
                             "unplaced-scaffold",
                             "pseudo-scaffold",
                             "fix-patch",
                             "novel-patch")
    sequence_role <- factor(ans[ , "SequenceRole"], levels=SequenceRole_levels)
    stopifnot(identical(is.na(sequence_role), is.na(ans[ , "SequenceRole"])))
    ans[ , "SequenceRole"] <- sequence_role
    oo <- order(as.integer(sequence_role))
    ans <- S4Vectors:::extract_data_frame_rows(ans, oo)
    Relationship_levels <- c("=", "<>")
    ans[ , "Relationship"] <- factor(ans[ , "Relationship"],
                                     levels=Relationship_levels)
    ans[ , "AssemblyUnit"] <- factor(ans[ , "AssemblyUnit"])
    UCSCStyleName <- ans[ , "UCSCStyleName"]
    if (!is.character(UCSCStyleName))
        UCSCStyleName <- as.character(UCSCStyleName)
    na_idx <- which(UCSCStyleName %in% "na")
    UCSCStyleName[na_idx] <- NA_character_
    ans[ , "UCSCStyleName"] <- UCSCStyleName
    circular <- make_circ_flags_from_circ_seqs(ans[ , "SequenceName"],
                                               circ_seqs=circ_seqs)
    stopifnot(all(ans[which(circular), "SequenceRole"] %in%
                  "assembled-molecule"))
    ans$circular <- circular
    ans
}

.NCBI_cached_chrom_info <- new.env(parent=emptyenv())

.get_NCBI_chrom_info_from_accession <- function(accession, circ_seqs=NULL,
    assembled.molecules.only=FALSE,
    assembly.units=NULL,
    recache=FALSE)
{
    ans <- .NCBI_cached_chrom_info[[accession]]
    if (is.null(ans) || recache) {
        assembly_report <- fetch_assembly_report(accession)
        ans <- .format_assembly_report(assembly_report, circ_seqs=circ_seqs)
        .NCBI_cached_chrom_info[[accession]] <- ans
    }
    if (assembled.molecules.only) {
        keep_idx <- which(ans[ , "SequenceRole"] %in% "assembled-molecule")
        ans <- S4Vectors:::extract_data_frame_rows(ans, keep_idx)
    }
    if (!is.null(assembly.units)) {
        lower_case_units <- tolower(assembly.units)
        ans_assembly_units <- ans[ , "AssemblyUnit"]
        ## "Primary Assembly" and "non-nuclear" are **always** considered valid
        ## Assembly Units regardless of whether they appear in the AssemblyUnit
        ## column or not.
        valid_assembly_units <- tolower(c("Primary Assembly",
                                          "non-nuclear",
                                          levels(ans_assembly_units)))
        bad_idx <- which(!(lower_case_units %in% valid_assembly_units))
        if (length(bad_idx) != 0L) {
            in1string <- paste0(assembly.units[bad_idx], collapse=", ")
            stop(wmsg("invalid Assembly Unit(s): ", in1string))
        }
        keep_idx <- which(tolower(ans_assembly_units) %in% lower_case_units)
        ans <- S4Vectors:::extract_data_frame_rows(ans, keep_idx)
    }
    ans
}

### Return an 8-column data.frame with columns "SequenceName" (character),
### "SequenceRole" (factor),  "GenBankAccn" (character), "Relationship"
### (factor), "RefSeqAccn" (character), "AssemblyUnit" (factor),
### "SequenceLength" (integer), "UCSCStyleName" (character),
### and "circular" (logical).
getChromInfoFromNCBI <- function(genome,
    assembled.molecules.only=FALSE,
    assembly.units=NULL,
    recache=FALSE,
    as.Seqinfo=FALSE)
{
    if (!isSingleString(genome))
        stop(wmsg("'genome' must be a single string"))
    if (!isTRUEorFALSE(assembled.molecules.only))
        stop(wmsg("'assembled.molecules.only' must be TRUE or FALSE"))
    if (!is.null(assembly.units)) {
        if (!is.character(assembly.units))
            stop(wmsg("'assembly.units' must be NULL or a character vector"))
        if (anyNA(assembly.units) ||
            !all(nzchar(assembly.units)) ||
            anyDuplicated(assembly.units))
            stop(wmsg("'assembly.units' cannot contain NAs, empty strings, ",
                      "or duplicates"))
    }
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(as.Seqinfo))
        stop(wmsg("'as.Seqinfo' must be TRUE or FALSE"))

    ## First see if the user supplied the assembly accession of a registered
    ## genome assembly instead of its name.
    assembly <- lookup_NCBI_accession2assembly(genome)
    if (!is.null(assembly)) {
        ## Yes s/he did.
        accession <- genome
        circ_seqs <- assembly$circ_seqs
    } else {
        ## No s/he didn't.
        ## Now see if s/he supplied the name of a registered genome assembly.
        accession <- lookup_NCBI_genome2accession(genome)
        if (!is.null(accession)) {
            ## Yes s/he did.
            if (length(accession) > 1L) {
                in1string <- paste0(accession, collapse=", ")
                warning(wmsg("Genome ", genome, " is mapped to more ",
                             "than one assembly (", in1string, "). ",
                             "The first one was selected."))
                accession <- accession[[1L]]
            }
            assembly <- lookup_NCBI_accession2assembly(accession)
            circ_seqs <- assembly$circ_seqs
        } else {
            ## No s/he didn't.
            ## So now we just assume that 'genome' is an assembly accession
            ## (an unregistered one).
            accession <- genome
            circ_seqs <- NULL
        }
    }
    ans <- .get_NCBI_chrom_info_from_accession(accession,
                circ_seqs=circ_seqs,
                assembled.molecules.only=assembled.molecules.only,
                assembly.units=assembly.units,
                recache=recache)
    if (!as.Seqinfo)
        return(ans)
    Seqinfo(seqnames=ans[ , "SequenceName"],
            seqlengths=ans[ , "SequenceLength"],
            isCircular=ans[ , "circular"],
            genome=genome)
}

