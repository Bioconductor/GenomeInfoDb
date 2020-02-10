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
                      "names not found in the Ensembl dataset"))
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
    special_mappings=NULL)
{
    NCBI_chrom_info <- getChromInfoFromNCBI(assembly_accession)

    Ensembl_seqlevels <- Ensembl_chrom_info[ , "name"]
    Ensembl_synonyms  <- Ensembl_chrom_info[ , "synonyms"]
    NCBI_seqlevels    <- NCBI_chrom_info[ , "SequenceName"]
    NCBI_GenBankAccn  <- NCBI_chrom_info[ , "GenBankAccn"]
    NCBI_RefSeqAccn   <- NCBI_chrom_info[ , "RefSeqAccn"]
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

.normarg_coord.systems <- function(coord.systems)
{
    if (is.null(coord.systems))
        return(NULL)
    if (!is.character(coord.systems))
        stop(wmsg("'coord.systems' must be a character vector or NULL"))
    stop_if_not_primary_key(coord.systems, "'coord.systems'")
    ## Sorting pushes normalization of the 'coord.systems' vector even
    ## further. This is important because we will encode the normalized vector
    ## in the caching key.
    sort(coord.systems)
}

.format_Ensembl_chrom_info <- function(seq_regions, circ_seqs=NULL)
{
    drop_columns <- c("seq_region_id", "coord_system_id",
                      "coord_system.species_id", "coord_system.version",
                      "coord_system.rank", "coord_system.attrib")
    ans <- drop_cols(seq_regions, drop_columns)
    ans <- rename_cols(ans, "coord_system.name", "coord_system")

    ## Column "coord_system".
    ans[ , "coord_system"] <- factor(ans[ , "coord_system"])

    ## Add column "circular".
    circular <- make_circ_flags_from_circ_seqs(ans[ , "name"],
                                               circ_seqs=circ_seqs)
    ans$circular <- circular

    ans
}

.ENSEMBL_cached_chrom_info <- new.env(parent=emptyenv())

.make_caching_key <- function(core_url, coord.systems)
    paste(c(core_url, coord.systems), collapse=":")

.get_chrom_info_for_unregistered_Ensembl_dataset <- function(core_url,
    assembled.molecules.only=FALSE,
    coord.systems=NULL,
    include.non_ref.sequences=FALSE,
    include.contigs=FALSE,
    include.clones=FALSE,
    recache=FALSE)
{
    coord.systems <- .normarg_coord.systems(coord.systems)
    caching_key <- .make_caching_key(core_url, coord.systems)
    ans <- .ENSEMBL_cached_chrom_info[[caching_key]]
    if (is.null(ans) || recache) {
        seq_regions <- fetch_seq_regions_from_Ensembl_ftp(core_url,
                                         coord_system_names=coord.systems,
                                         add.toplevel.col=TRUE,
                                         add.non_ref.col=TRUE)
        ans <- .format_Ensembl_chrom_info(seq_regions)
        .ENSEMBL_cached_chrom_info[[caching_key]] <- ans
    }
    if (assembled.molecules.only) {
        ## FIXME: This is broken on some datasets e.g. btaurus_gene_ensembl
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
    ans
}

getChromInfoFromEnsembl <- function(dataset,
    release=NA, use.grch37=FALSE, kingdom=NA,
    assembled.molecules.only=FALSE,
    coord.systems=NULL,
    include.non_ref.sequences=FALSE,
    include.contigs=FALSE,
    include.clones=FALSE,
    recache=FALSE,
    as.Seqinfo=FALSE)
{
    ## TODO: get_url_to_Ensembl_ftp_mysql_core() is not users proof (i.e.
    ## it does not check user input). Either make it users proof or wrap it
    ## into a users proof wrapper.
    core_url <- get_url_to_Ensembl_ftp_mysql_core(dataset,
                                                  release=release,
                                                  use.grch37=use.grch37,
                                                  kingdom=kingdom)
    if (!isTRUEorFALSE(assembled.molecules.only))
        stop(wmsg("'assembled.molecules.only' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(include.non_ref.sequences))
        stop(wmsg("'include.non_ref.sequences' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(include.contigs))
        stop(wmsg("'include.contigs' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(include.clones))
        stop(wmsg("'include.clones' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(as.Seqinfo))
        stop(wmsg("'as.Seqinfo' must be TRUE or FALSE"))

    ans <- .get_chrom_info_for_unregistered_Ensembl_dataset(core_url,
                assembled.molecules.only=assembled.molecules.only,
                coord.systems=coord.systems,
                include.non_ref.sequences=include.non_ref.sequences,
                include.contigs=include.contigs,
                include.clones=include.clones,
                recache=recache)

    if (!as.Seqinfo)
        return(ans)
    Seqinfo(seqnames=ans[ , "name"],
            seqlengths=ans[ , "length"],
            isCircular=ans[ , "circular"])
}

