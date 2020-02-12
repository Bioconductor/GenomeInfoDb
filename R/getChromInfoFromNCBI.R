### =========================================================================
### getChromInfoFromNCBI()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load and access db of registered NCBI assemblies
###

.NCBI_assembly2accession <- new.env(parent=emptyenv())
.NCBI_accession2assembly_info <- new.env(parent=emptyenv())

.load_registered_NCBI_assembly <- function(file_path)
{
    filename <- basename(file_path)
    if (substr(filename, nchar(filename)-1L, nchar(filename)) != ".R")
        stop(wmsg("name of assembly registration file '", filename, "' must ",
                  "have extension .R"))
    if (grepl(" ", filename, fixed=TRUE))
        stop(wmsg("name of assembly registration file '", filename, "' must ",
                  "not contain spaces"))

    ## Placeholders. Will actually get defined when we source the
    ## assembly files.
    ORGANISM <- NULL    # Expected to be a single string.
    ASSEMBLIES <- NULL  # Expected to be a list with one list element per
                        # assembly.
    source(file_path, local=TRUE)

    stop_if <- function(notok, ...) {
        if (notok)
            stop("Error in NCBI assembly registration file '", filename,
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

    required_fields <- c("assembly", "date", "assembly_accession", "circ_seqs")
    for (i in seq_along(ASSEMBLIES)) {
        assembly_info <- ASSEMBLIES[[i]]
        label <- sprintf("'ASSEMBLIES[[%d]]'", i)

        stop_if(!is.list(assembly_info), label, " must be a named list")
        assembly_fields <- names(assembly_info)
        stop_if(!is.character(assembly_fields), label, " must be a named list")
        stop_if(!is_primary_key(assembly_fields),
                "the names on ", label, " must ",
                "not contain NAs, empty strings, or duplicates")
        stop_if(!all(required_fields %in% names(assembly_info)),
                label, " must have fields: ",
                paste(paste0("\"", required_fields, "\""), collapse=", "))

        ## Check "assembly" field (required).
        assembly <- assembly_info$assembly
        stop_if(!isSingleString(assembly) || assembly == "",
                "\"assembly\" field in ", label, " must ",
                "be a single non-empty string")

        ## Check "date" field (required).
        date <- assembly_info$date
        stop_if(!isSingleString(date) || date == "",
                "\"date\" field in ", label, " must ",
                "be a single non-empty string")

        ## Check "accession" field (required).
        accession <- assembly_info$assembly_accession
        stop_if(!isSingleString(accession) || accession == "",
                "\"accession\" field in ", label, " must ",
                "be a single non-empty string")
        stop_if(!is.null(.NCBI_accession2assembly_info[[accession]]),
                "assembly accession ", accession, " used more than once")

        ## Check "circ_seqs" field (required).
        circ_seqs <- assembly_info$circ_seqs
        stop_if(!is.character(circ_seqs),
                "\"circ_seqs\" field in ", label, " must ",
                "be a character vector")
        stop_if(!is_primary_key(circ_seqs),
                "\"circ_seqs\" field in ", label, " must ",
                "not contain NAs, empty strings, or duplicates")

        ## Check optional fields.

        extra_info <- assembly_info$extra_info
        if (!is.null(extra_info)) {
            stop_if(!is.character(extra_info),
                "\"extra_info\" field in ", label, " must ",
                "be a named character vector")
            stop_if(anyNA(extra_info) || !all(nzchar(extra_info)),
                "\"extra_info\" field in ", label, " must ",
                "not contain NAs or empty strings")
            tags <- names(extra_info)
            stop_if(!is.character(tags),
                "\"extra_info\" field in ", label, " must ",
                "be a named character vector")
            stop_if(anyNA(tags) || anyDuplicated(tags),
                "the names on \"extra_info\" field ", label, " must ",
                "not contain NAs or duplicates")
        }

        special_mappings <- assembly_info$NCBI2Ensembl_special_mappings
        if (!is.null(special_mappings)) {
            stop_if(!is.character(special_mappings),
                "\"NCBI2Ensembl_special_mappings\" field in ", label, " must ",
                "be a named character vector")
            stop_if(is.null(names(special_mappings)),
                "\"NCBI2Ensembl_special_mappings\" field in ", label, " must ",
                "be a named character vector")
        }

        assembly <- tolower(assembly)
        .NCBI_assembly2accession[[assembly]] <-
            c(.NCBI_assembly2accession[[assembly]], accession)
        assembly_info$organism <- ORGANISM
        assembly_info$rank_within_organism <- i
        .NCBI_accession2assembly_info[[accession]] <- assembly_info
    }
}

.load_registered_NCBI_assemblies <- function()
{
    dir_path <- system.file("registered", "NCBI_assemblies",
                             package="GenomeInfoDb")
    file_paths <- list.files(dir_path, pattern="\\.R$", full.names=TRUE)
    for (file_path in file_paths)
        .load_registered_NCBI_assembly(file_path)
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

registered_NCBI_assemblies <- function()
{
    if (length(.NCBI_accession2assembly_info) == 0L)
        .load_registered_NCBI_assemblies()
    assemblies <- unname(as.list(.NCBI_accession2assembly_info, all.names=TRUE))

    colnames <- c("organism", "rank_within_organism", "assembly", "date",
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

.lookup_NCBI_assembly2accession <- function(assembly)
{
    stopifnot(isSingleString(assembly))
    if (length(.NCBI_assembly2accession) == 0L)
        .load_registered_NCBI_assemblies()
    .NCBI_assembly2accession[[tolower(assembly)]]
}

find_NCBI_assembly_info_for_accession <- function(accession)
{
    stopifnot(isSingleString(accession))
    if (length(.NCBI_accession2assembly_info) == 0L)
        .load_registered_NCBI_assemblies()
    .NCBI_accession2assembly_info[[accession]]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getChromInfoFromNCBI()
###

.format_NCBI_chrom_info <- function(assembly_report, circ_seqs=NULL)
{
    ans <- drop_cols(assembly_report, "AssignedMoleculeLocationOrType")

    ## Column "SequenceName".
    sequence_name <- as.character(ans[ , "SequenceName"])
    stopifnot(is_primary_key(sequence_name))
    ans[ , "SequenceName"] <- sequence_name

    ## Column "SequenceRole".
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

    ## Re-order the rows based on SequenceRole.
    oo <- order(as.integer(sequence_role))
    ans <- S4Vectors:::extract_data_frame_rows(ans, oo)

    ## Column "AssignedMolecule".
    idx <- which(ans[ , "SequenceRole"] %in% "assembled-molecule")
    assembled_molecules <- ans[idx , "SequenceName"]
    ans[ , "AssignedMolecule"] <- factor(ans[ , "AssignedMolecule"],
                                         levels=assembled_molecules)

    ## Column "Relationship".
    Relationship_levels <- c("=", "<>")
    ans[ , "Relationship"] <- factor(ans[ , "Relationship"],
                                     levels=Relationship_levels)

    ## Column "AssemblyUnit".
    ans[ , "AssemblyUnit"] <- factor(ans[ , "AssemblyUnit"])

    ## Column "UCSCStyleName".
    UCSCStyleName <- ans[ , "UCSCStyleName"]
    if (!is.character(UCSCStyleName))
        UCSCStyleName <- as.character(UCSCStyleName)
    na_idx <- which(UCSCStyleName %in% "na")
    UCSCStyleName[na_idx] <- NA_character_
    ans[ , "UCSCStyleName"] <- UCSCStyleName

    ## Add column "circular".
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
        ans <- .format_NCBI_chrom_info(assembly_report, circ_seqs=circ_seqs)
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

### Return a 10-column data.frame with columns "SequenceName" (character),
### "SequenceRole" (factor),  "AssignedMolecule" (factor), "GenBankAccn"
### (character), "Relationship" (factor), "RefSeqAccn" (character),
### "AssemblyUnit" (factor), "SequenceLength" (integer), "UCSCStyleName"
### (character), and "circular" (logical).
getChromInfoFromNCBI <- function(assembly,
    assembled.molecules.only=FALSE,
    assembly.units=NULL,
    recache=FALSE,
    as.Seqinfo=FALSE)
{
    if (!isSingleString(assembly))
        stop(wmsg("'assembly' must be a single string"))
    if (!isTRUEorFALSE(assembled.molecules.only))
        stop(wmsg("'assembled.molecules.only' must be TRUE or FALSE"))
    if (!is.null(assembly.units)) {
        if (!is.character(assembly.units))
            stop(wmsg("'assembly.units' must be a character vector or NULL"))
        stop_if_not_primary_key(assembly.units, "'assembly.units'")
    }
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(as.Seqinfo))
        stop(wmsg("'as.Seqinfo' must be TRUE or FALSE"))

    ## First see if the user supplied the accession of a registered assembly
    ## instead of the name of the assembly.
    NCBI_assembly_info <- find_NCBI_assembly_info_for_accession(assembly)
    if (!is.null(NCBI_assembly_info)) {
        ## Yes s/he did.
        accession <- assembly
        circ_seqs <- NCBI_assembly_info$circ_seqs
    } else {
        ## No s/he didn't.
        ## Now see if s/he supplied the name of a registered assembly.
        accession <- .lookup_NCBI_assembly2accession(assembly)
        if (!is.null(accession)) {
            ## Yes s/he did.
            if (length(accession) > 1L) {
                in1string <- paste0(accession, collapse=", ")
                warning(wmsg("Assembly ", assembly, " is mapped to more ",
                             "than one assembly (", in1string, "). ",
                             "The first one was selected."))
                accession <- accession[[1L]]
            }
            NCBI_assembly_info <-
                find_NCBI_assembly_info_for_accession(accession)
            circ_seqs <- NCBI_assembly_info$circ_seqs
        } else {
            ## No s/he didn't.
            ## So now we just assume that 'assembly' is an assembly accession
            ## (an unregistered one).
            accession <- assembly
            circ_seqs <- NULL
        }
    }

    ans <- .get_NCBI_chrom_info_from_accession(accession,
                circ_seqs=circ_seqs,
                assembled.molecules.only=assembled.molecules.only,
                assembly.units=assembly.units,
                recache=recache)

    if (!as.Seqinfo) {
        attr(ans, "NCBI_assembly_info") <- NCBI_assembly_info
        return(ans)
    }
    if (is.null(NCBI_assembly_info)) {
        ans_genome <- assembly
    } else {
        ans_genome <- NCBI_assembly_info$assembly
    }
    Seqinfo(seqnames=ans[ , "SequenceName"],
            seqlengths=ans[ , "SequenceLength"],
            isCircular=ans[ , "circular"],
            genome=ans_genome)
}

