### =========================================================================
### getChromInfoFromEnsembl()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .add_NCBI_cols_to_Ensembl_chrom_info()
###

.match_Ensembl_synonmys_to_NCBI_col <- function(Ensembl_seqlevels,
                                                Ensembl_synonyms,
                                                NCBI_accns,
                                                what)
{
    stopifnot(is(Ensembl_synonyms, "CompressedCharacterList"))
    Ensembl_unlisted_synonyms <- unlist(Ensembl_synonyms, use.names=FALSE)

    ## solid_match() will fail with a not-so-useful error message
    ## if 'Ensembl_unlisted_synonyms' contains duplicates. We check
    ## this early so we can display a slightly better error message.
    if (anyDuplicated(Ensembl_unlisted_synonyms))
        stop(wmsg("'unlist(Ensembl_synonyms)' contains duplicates"))

    m <- solid_match(unname(Ensembl_unlisted_synonyms), NCBI_accns,
                     x_what="Ensembl synonym",
                     table_what=what)

    m <- relist(m, Ensembl_synonyms)
    m <- m[!is.na(m)]
    m_lens <- lengths(m)
    ambig_idx <- which(m_lens > 1L)
    if (length(ambig_idx) != 0L) {
        in1string <- paste0(Ensembl_seqlevels[ambig_idx], collapse = ", ")
        stop(wmsg("Ensembl seqlevel(s) matched to more ",
                  "than 1 ", what, ": ", in1string))
    }
    as.integer(m)
}

.mk_progress_str <- function(L2R, NCBI_seqlevels)
{
    N <- length(NCBI_seqlevels)
    ## Using 'sum(!is.na(L2R))' would work but it is based on the assumption
    ## that 'L2R' will never contain duplicates. We don't want to make this
    ## assumption here so let's stay on the safe side.
    revmap_count <- sum(tabulate(L2R, nbins=length(NCBI_seqlevels)) != 0L)
    percent <- 100 * revmap_count / N
    sprintf("Total NCBI seqlevels mapped: %d/%d -- %.2f%%",
            revmap_count, N, percent)
}

### The workhorse behind .add_NCBI_cols_to_Ensembl_chrom_info().
### - All input vectors (except 'Ensembl_synonyms') must be character vectors.
### - 'Ensembl_synonyms' must be a list-like object parallel to
###   'Ensembl_seqlevels' where each list element is a character vector
###   containing all the synonyms for the corresponding Ensembl seqlevel.
### - All NCBI input vectors must have the same length.
### - Vectors 'Ensembl_seqlevels' and 'NCBI_seqlevels' must be "primary
###   keys" i.e. must not contain NAs, empty strings, or duplicates.
###   (No such assumption is made about the other input vectors.)
### Returns an integer vector parallel to 'Ensembl_seqlevels'.
.map_Ensembl_seqlevels_to_NCBI_seqlevels <- function(Ensembl_seqlevels,
                                                     Ensembl_synonyms,
                                                     NCBI_seqlevels,
                                                     NCBI_GenBankAccn,
                                                     NCBI_RefSeqAccn,
                                                     special_mappings=NULL,
                                                     verbose=FALSE)
{
    stopifnot(is.character(Ensembl_seqlevels),
              is_primary_key(Ensembl_seqlevels),
              is(Ensembl_synonyms, "list_OR_List"),
              is.character(NCBI_seqlevels),
              is_primary_key(NCBI_seqlevels),
              is.character(NCBI_GenBankAccn),
              is.character(NCBI_RefSeqAccn))
    if (!is(Ensembl_synonyms, "CompressedCharacterList"))
        Ensembl_synonyms <- as(Ensembl_synonyms, "CompressedCharacterList")
    stopifnot(length(Ensembl_synonyms) == length(Ensembl_seqlevels))

    prefix <- "Mapping Ensembl to NCBI: "
    L2R <- rep.int(NA_integer_, length(Ensembl_seqlevels))
    unmapped_idx <- seq_along(L2R)

    ## 1. Handle special mappings.
    if (!is.null(special_mappings)) {
        if (verbose)
            message(prefix, "Handle special mappings ... ", appendLF=FALSE)
        m1 <- match(names(special_mappings), Ensembl_seqlevels)
        if (anyNA(m1))
            stop(wmsg("'names(special_mappings)' contains sequence ",
                      "names not found in Ensembl core db"))
        m2 <- match(special_mappings, NCBI_seqlevels)
        if (anyNA(m2))
            stop(wmsg("'special_mappings' contains sequence ",
                      "names not found in the NCBI assembly"))
        L2R[m1] <- m2
        if (verbose)
            message("OK (", .mk_progress_str(L2R, NCBI_seqlevels), ")")

        unmapped_idx <- which(is.na(L2R))
        if (length(unmapped_idx) == 0L)
            return(L2R)
    }

    ## 2. Assign based on exact matching (case insensitive) of
    ##    the seqlevels.
    if (verbose)
        message(prefix, "Ensembl_seqlevels -> NCBI_seqlevels ... ",
                appendLF=FALSE)
    m <- solid_match2(Ensembl_seqlevels, NCBI_seqlevels,
                      x_what="Ensembl seqlevel",
                      table_what="NCBI seqlevel")
    L2R[unmapped_idx] <- m[unmapped_idx]
    if (verbose)
        message("OK (", .mk_progress_str(L2R, NCBI_seqlevels), ")")

    unmapped_idx <- which(is.na(L2R))
    if (length(unmapped_idx) == 0L)
        return(L2R)

    ## 3. Match 'Ensembl_seqlevels' to 'NCBI_GenBankAccn'.
    if (verbose)
        message(prefix, "Ensembl_seqlevels -> NCBI_GenBankAccn ... ",
                appendLF=FALSE)
    m <- solid_match2(Ensembl_seqlevels, NCBI_GenBankAccn,
                      x_what="Ensembl seqlevel",
                      table_what="NCBI GenBankAccn value")
    L2R[unmapped_idx] <- m[unmapped_idx]
    if (verbose)
        message("OK (", .mk_progress_str(L2R, NCBI_seqlevels), ")")

    unmapped_idx <- which(is.na(L2R))
    if (length(unmapped_idx) == 0L)
        return(L2R)

    ## 4. Match 'Ensembl_seqlevels' to 'NCBI_RefSeqAccn'.
    if (verbose)
        message(prefix, "Ensembl_seqlevels -> NCBI_RefSeqAccn ... ",
                appendLF=FALSE)
    m <- solid_match2(Ensembl_seqlevels, NCBI_RefSeqAccn,
                      x_what="Ensembl seqlevel",
                      table_what="NCBI RefSeqAccn value")
    L2R[unmapped_idx] <- m[unmapped_idx]
    if (verbose)
        message("OK (", .mk_progress_str(L2R, NCBI_seqlevels), ")")

    unmapped_idx <- which(is.na(L2R))
    if (length(unmapped_idx) == 0L)
        return(L2R)

    ## 5. Match 'Ensembl_synonyms' to 'NCBI_seqlevels'.
    if (verbose)
        message(prefix, "Ensembl_synonyms -> NCBI_seqlevels ... ",
                appendLF=FALSE)
    m <- .match_Ensembl_synonmys_to_NCBI_col(Ensembl_seqlevels,
                                             Ensembl_synonyms,
                                             NCBI_seqlevels,
                                             what="NCBI seqlevel")
    L2R[unmapped_idx] <- m[unmapped_idx]
    if (verbose)
        message("OK (", .mk_progress_str(L2R, NCBI_seqlevels), ")")

    unmapped_idx <- which(is.na(L2R))
    if (length(unmapped_idx) == 0L)
        return(L2R)

    ## 6. Match 'Ensembl_synonyms' to 'NCBI_GenBankAccn'.
    if (verbose)
        message(prefix, "Ensembl_synonyms -> NCBI_GenBankAccn ... ",
                appendLF=FALSE)
    m <- .match_Ensembl_synonmys_to_NCBI_col(Ensembl_seqlevels,
                                             Ensembl_synonyms,
                                             NCBI_GenBankAccn,
                                             what="NCBI GenBankAccn value")
    L2R[unmapped_idx] <- m[unmapped_idx]
    if (verbose)
        message("OK (", .mk_progress_str(L2R, NCBI_seqlevels), ")")

    unmapped_idx <- which(is.na(L2R))
    if (length(unmapped_idx) == 0L)
        return(L2R)

    ## 7. Match 'Ensembl_synonyms' to 'NCBI_RefSeqAccn'.
    if (verbose)
        message(prefix, "Ensembl_synonyms -> NCBI_RefSeqAccn ... ",
                appendLF=FALSE)
    m <- .match_Ensembl_synonmys_to_NCBI_col(Ensembl_seqlevels,
                                             Ensembl_synonyms,
                                             NCBI_RefSeqAccn,
                                             what="NCBI RefSeqAccn value")
    L2R[unmapped_idx] <- m[unmapped_idx]
    if (verbose)
        message("OK (", .mk_progress_str(L2R, NCBI_seqlevels), ")")

    L2R
}

.add_NCBI_cols_to_Ensembl_chrom_info <- function(
    Ensembl_chrom_info, assembly_accession,
    NCBI2Ensembl_special_mappings=NULL)
{
    NCBI_chrom_info <- getChromInfoFromNCBI(assembly_accession)

    Ensembl_seqlevels <- Ensembl_chrom_info[ , "name"]
    Ensembl_synonyms  <- Ensembl_chrom_info[ , "synonyms"]
    NCBI_seqlevels    <- NCBI_chrom_info[ , "SequenceName"]
    NCBI_GenBankAccn  <- NCBI_chrom_info[ , "GenBankAccn"]
    NCBI_RefSeqAccn   <- NCBI_chrom_info[ , "RefSeqAccn"]

    ## Reverse the mappings.
    if (is.null(NCBI2Ensembl_special_mappings)) {
        special_mappings <- NULL
    } else {
        special_mappings <- names(NCBI2Ensembl_special_mappings)
        names(special_mappings) <- NCBI2Ensembl_special_mappings
    }

    L2R <- .map_Ensembl_seqlevels_to_NCBI_seqlevels(
                                          Ensembl_seqlevels,
                                          Ensembl_synonyms,
                                          NCBI_seqlevels,
                                          NCBI_GenBankAccn,
                                          NCBI_RefSeqAccn,
                                          special_mappings=special_mappings)

    NCBI_chrom_info <- S4Vectors:::extract_data_frame_rows(NCBI_chrom_info, L2R)

    compare_idx <- which(!is.na(NCBI_chrom_info[ , "SequenceLength"]))
    stopifnot(identical(Ensembl_chrom_info[compare_idx, "length"],
                        NCBI_chrom_info[compare_idx, "SequenceLength"]))

    drop_colnames <- c("SequenceLength", "UCSCStyleName", "circular")
    cbind(Ensembl_chrom_info, drop_cols(NCBI_chrom_info, drop_colnames))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getChromInfoFromEnsembl()
###

.format_Ensembl_chrom_info <- function(seq_regions, circ_seqs=NULL)
{
    drop_columns <- c("seq_region_id", "coord_system_id",
                      "coord_system.species_id", "coord_system.version",
                      "coord_system.rank", "coord_system.attrib")
    ans <- drop_cols(seq_regions, drop_columns)

    ## Column "coord_system".
    ans <- rename_cols(ans, "coord_system.name", "coord_system")
    ans[ , "coord_system"] <- factor(ans[ , "coord_system"])

    ## Add column "circular".
    circular <- make_circ_flags_from_circ_seqs(ans[ , "name"],
                                               circ_seqs=circ_seqs)
    ans$circular <- circular

    ans
}

### Return NULL or a named list with at least component "assembly_accession"
### and optionally component "NCBI2Ensembl_special_mappings" (plus possibly
### other components that are will be ignored by the caller).
.find_NCBI_assembly_info_for_Ensembl_species_info <- function(species_info)
{
    check_species_info(species_info)
    species <- species_info$species
    Ensembl_release <- species_info$Ensembl_release

    assembly_accession <- species_info$assembly_accession
    if (is.null(assembly_accession)) {
        warning(wmsg("'map.NCBI' got ignored for species \"",
                     species, "\" (no assembly accession can be found ",
                     "for it in Ensembl release ", Ensembl_release, ")"))
        return(NULL)
    }

    ## Is this Ensembl species associated with a registered NCBI assembly?

    NCBI_assembly_info <-
        find_NCBI_assembly_info_for_accession(assembly_accession)
    if (!is.null(NCBI_assembly_info))
        return(NCBI_assembly_info)  # yes!

    fallback_assembly_info <- list(assembly_accession=assembly_accession)

    Ensembl_assembly <- species_info$assembly
    if (is.null(Ensembl_assembly))
        return(fallback_assembly_info)  # no

    NCBI_assemblies <- registered_NCBI_assemblies()
    idx <- which(NCBI_assemblies[ , "assembly"] %in% Ensembl_assembly)
    if (length(idx) == 0L)
        return(fallback_assembly_info)  # no

    if (length(idx) == 1L) {        # yes!
        accession <- NCBI_assemblies[idx , "assembly_accession"]
        NCBI_assembly_info <- find_NCBI_assembly_info_for_accession(accession)
        return(NCBI_assembly_info)
    }
    ## Yes, but ambiguously! We stick to the accession provided by Ensembl
    ## with a warning.
    in1string <- paste(NCBI_assemblies[idx , "assembly_accession"],
                       collapse=", ")
    warning(wmsg(
        "Assembly accession found for species \"", species, "\" in ",
        "Ensembl release ", Ensembl_release, " is ", assembly_accession, ". ",
        "This is not associated with **any** of the NCBI assemblies ",
        "registered in GenomeInfoDb (see registered_NCBI_assemblies()). ",
        "However, using the assembly name (", Ensembl_assembly, ") ",
        "provided by Ensembl in release ", Ensembl_release, ", ",
        "it can be associated with **more than one** registered ",
        "NCBI assembly (", in1string, "). Because this association is ",
        "ambiguous, we ignored it and associated \"", species, "\" ",
        "with the assembly accession provided by Ensembl."))
    fallback_assembly_info
}

.ENSEMBL_cached_chrom_info <- new.env(parent=emptyenv())

.get_chrom_info_from_Ensembl_FTP <- function(core_db_url,
    assembled.molecules.only=FALSE,
    include.non_ref.sequences=FALSE,
    include.contigs=FALSE,
    include.clones=FALSE,
    species_info=NULL,
    recache=FALSE)
{
    ans <- .ENSEMBL_cached_chrom_info[[core_db_url]]
    if (is.null(ans) || recache) {
        seq_regions <- fetch_seq_regions_from_Ensembl_FTP(core_db_url,
                                         add.toplevel.col=TRUE,
                                         add.non_ref.col=TRUE)
        ans <- .format_Ensembl_chrom_info(seq_regions)
        .ENSEMBL_cached_chrom_info[[core_db_url]] <- ans
    }

    ## Add NCBI cols.
    if (!is.null(species_info)) {
        NCBI_assembly_info <-
            .find_NCBI_assembly_info_for_Ensembl_species_info(species_info)
        if (!is.null(NCBI_assembly_info)) {
            ans <- .add_NCBI_cols_to_Ensembl_chrom_info(
                             ans, NCBI_assembly_info$assembly_accession,
                             NCBI_assembly_info$NCBI2Ensembl_special_mappings)
        }
    }

    if (assembled.molecules.only) {
        ## FIXME: This is broken for some core dbs e.g. bos_taurus_core_99_12
        ## where coord_system is not set to "chromosome" for chromosomes.
        keep_idx <- which(ans[ , "coord_system"] %in% "chromosome")
        ans <- S4Vectors:::extract_data_frame_rows(ans, keep_idx)
    }
    if (!include.non_ref.sequences) {
        keep_idx <- which(!ans[ , "non_ref"])
        ans <- S4Vectors:::extract_data_frame_rows(ans, keep_idx)
    }
    if (!include.contigs) {
        keep_idx <- which(!(ans[ , "coord_system"] %in% "contig"))
        ans <- S4Vectors:::extract_data_frame_rows(ans, keep_idx)
    }
    if (!include.clones) {
        keep_idx <- which(!(ans[ , "coord_system"] %in% "clone"))
        ans <- S4Vectors:::extract_data_frame_rows(ans, keep_idx)
    }

    if (!(is.null(species_info) || is.null(NCBI_assembly_info)))
        attr(ans, "NCBI_assembly_info") <- NCBI_assembly_info
    ans
}

getChromInfoFromEnsembl <- function(species,
    release=NA, division=NA, use.grch37=FALSE,
    assembled.molecules.only=FALSE,
    include.non_ref.sequences=FALSE,
    include.contigs=FALSE,
    include.clones=FALSE,
    map.NCBI=FALSE,
    recache=FALSE,
    as.Seqinfo=FALSE)
{
    core_db_url <- get_Ensembl_FTP_core_db_url(species, release=release,
                                               division=division,
                                               use.grch37=use.grch37)
    species_info <- attr(core_db_url, "species_info")
    check_species_info(species_info)

    if (!isTRUEorFALSE(assembled.molecules.only))
        stop(wmsg("'assembled.molecules.only' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(include.non_ref.sequences))
        stop(wmsg("'include.non_ref.sequences' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(include.contigs))
        stop(wmsg("'include.contigs' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(include.clones))
        stop(wmsg("'include.clones' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(map.NCBI))
        stop(wmsg("'map.NCBI' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(as.Seqinfo))
        stop(wmsg("'as.Seqinfo' must be TRUE or FALSE"))

    ans <- .get_chrom_info_from_Ensembl_FTP(core_db_url,
                assembled.molecules.only=assembled.molecules.only,
                include.non_ref.sequences=include.non_ref.sequences,
                include.contigs=include.contigs,
                include.clones=include.clones,
                species_info=if (map.NCBI) species_info else NULL,
                recache=recache)

    if (!as.Seqinfo) {
        attr(ans, "species_info") <- species_info
        return(ans)
    }
    ans_genome <- NA_character_
    assembly <- species_info$assembly
    if (!is.null(assembly))
        ans_genome <- assembly
    Seqinfo(seqnames=ans[ , "name"],
            seqlengths=ans[ , "length"],
            isCircular=ans[ , "circular"],
            genome=ans_genome)
}

