### =========================================================================
### getChromInfoFromEnsembl()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .add_NCBI_cols_to_Ensembl_chrom_info()
###

.mk_progress_str <- function(L2R, NCBI_seqlevels)
{
    N <- length(NCBI_seqlevels)
    ## Using 'sum(!is.na(L2R))' would work but it is based on the assumption
    ## that 'L2R' will never contain duplicates. We don't want to make this
    ## assumption here so let's stay on the safe side.
    revmap_count <- sum(tabulate(L2R, nbins=length(NCBI_seqlevels)) != 0L)
    percent <- 100 * revmap_count / N
    sprintf("Total NCBI seqlevels mapped: %d/%d -- %.1f%%",
            revmap_count, N, percent)
}

### The workhorse behind .add_NCBI_cols_to_Ensembl_chrom_info().
### - All input vectors must be character vectors.
### - All Ensembl input vectors must have the same length.
### - All NCBI input vectors must have the same length.
### - Vectors 'Ensembl_seqlevels' and 'NCBI_seqlevels' must be "primary
###   keys" i.e. must not contain NAs, empty strings, or duplicates.
### No assumptions are made about the other input vectors.
### Returns an integer vector parallel to 'Ensembl_seqlevels'.
.map_Ensembl_seqlevels_to_NCBI_seqlevels <- function(Ensembl_seqlevels,
                                                     Ensembl_INSDC,
                                                     Ensembl_RefSeq,
                                                     NCBI_seqlevels,
                                                     NCBI_GenBankAccn,
                                                     NCBI_RefSeqAccn,
                                                     special_mappings=NULL,
                                                     verbose=FALSE)
{
    stopifnot(is.character(Ensembl_seqlevels),
              is_primary_key(Ensembl_seqlevels),
              is.character(Ensembl_INSDC),
              is.character(Ensembl_RefSeq),
              is.character(NCBI_seqlevels),
              is_primary_key(NCBI_seqlevels),
              is.character(NCBI_GenBankAccn),
              is.character(NCBI_RefSeqAccn))

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

    ## From now on matching will be case insensitive.
    lower_with_names <- function(x) setNames(tolower(x), x)
    ens_seqlevels <- lower_with_names(Ensembl_seqlevels)
    ens_insdc <- lower_with_names(Ensembl_INSDC)
    ens_refseq <- lower_with_names(Ensembl_RefSeq)
    ncbi_seqlevels <- lower_with_names(NCBI_seqlevels)
    ncbi_genbankaccn <- lower_with_names(NCBI_GenBankAccn)
    ncbi_refseqaccn <- lower_with_names(NCBI_RefSeqAccn)

    ## 2. Assign based on exact matching (case insensitive) of
    ##    the seqlevels.
    if (verbose)
        message(prefix, "Ensembl seqlevels -> NCBI seqlevels ... ",
                appendLF=FALSE)
    m <- solid_match2(ens_seqlevels, ncbi_seqlevels,
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
        message(prefix, "Ensembl seqlevels -> NCBI GenBankAccn ... ",
                appendLF=FALSE)
    m <- solid_match2(ens_seqlevels, ncbi_genbankaccn,
                      x_what="Ensembl seqlevel",
                      table_what="NCBI GenBankAccn value")
    L2R[unmapped_idx] <- m[unmapped_idx]
    if (verbose)
        message("OK (", .mk_progress_str(L2R, NCBI_seqlevels), ")")

    unmapped_idx <- which(is.na(L2R))
    if (length(unmapped_idx) == 0L)
        return(L2R)

    ## 4. Match 'Ensembl_INSDC' to 'NCBI_GenBankAccn'.
    if (verbose)
        message(prefix, "Ensembl_INSDC -> NCBI GenBankAccn ... ",
                appendLF=FALSE)
    m <- solid_match2(ens_insdc, ncbi_genbankaccn,
                      x_what="INSDC synonym",
                      table_what="NCBI GenBankAccn value")
    L2R[unmapped_idx] <- m[unmapped_idx]
    if (verbose)
        message("OK (", .mk_progress_str(L2R, NCBI_seqlevels), ")")

    unmapped_idx <- which(is.na(L2R))
    if (length(unmapped_idx) == 0L)
        return(L2R)

    ## 5. Match 'Ensembl_RefSeq' to 'NCBI_RefSeqAccn'.
    if (verbose)
        message(prefix, "Ensembl_RefSeq -> NCBI RefSeqAccn ... ",
                appendLF=FALSE)
    m <- solid_match2(ens_refseq, ncbi_refseqaccn,
                      x_what="RefSeq synonym",
                      table_what="NCBI RefSeqAccn value")
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
    Ensembl_INSDC     <- Ensembl_chrom_info[ , "INSDC"]
    Ensembl_RefSeq    <- Ensembl_chrom_info[ , "RefSeq"]
    NCBI_seqlevels    <- NCBI_chrom_info[ , "SequenceName"]
    NCBI_GenBankAccn  <- NCBI_chrom_info[ , "GenBankAccn"]
    NCBI_RefSeqAccn   <- NCBI_chrom_info[ , "RefSeqAccn"]
    L2R <- .map_Ensembl_seqlevels_to_NCBI_seqlevels(
                                          Ensembl_seqlevels,
                                          Ensembl_INSDC,
                                          Ensembl_RefSeq,
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

.normarg_coord.system.names <- function(coord.system.names)
{
    if (is.null(coord.system.names))
        return(NULL)
    if (!is.character(coord.system.names))
        stop(wmsg("'coord.system.names' must be a character vector or NULL"))
    stop_if_not_primary_key(coord.system.names, "'coord.system.names'")
    ## Sorting pushes normalization of the 'coord.system.names' vector even
    ## further. This is important because we will encode the normalized vector
    ## in the caching key.
    sort(coord.system.names)
}

.format_Ensembl_chrom_info <- function(seq_regions, circ_seqs=NULL)
{
    drop_columns <- c("seq_region_id", "coord_system_id",
                      "coord_system.species_id", "coord_system.version",
                      "coord_system.rank", "coord_system.attrib")
    ans <- drop_cols(seq_regions, drop_columns)

    circular <- make_circ_flags_from_circ_seqs(ans[ , "name"],
                                               circ_seqs=circ_seqs)
    ans$circular <- circular

    ans
}

.ENSEMBL_cached_chrom_info <- new.env(parent=emptyenv())

.make_caching_key <- function(core_url, coord.system.names)
    paste(c(core_url, coord.system.names), collapse=":")

.get_chrom_info_for_unregistered_Ensembl_dataset <- function(core_url,
    assembled.molecules.only=FALSE,
    coord.system.names=c("chromosome", "scaffold"),
    include.no_ref.sequences=FALSE,
    recache=FALSE)
{
    coord.system.names <- .normarg_coord.system.names(coord.system.names)
    caching_key <- .make_caching_key(core_url, coord.system.names)
    ans <- .ENSEMBL_cached_chrom_info[[caching_key]]
    if (is.null(ans) || recache) {
        seq_regions <- fetch_seq_regions_from_Ensembl_ftp(core_url,
                                         coord_system_names=coord.system.names,
                                         add.toplevel.col=TRUE,
                                         add.non_ref.col=TRUE,
                                         add.INSDC.col=TRUE,
                                         add.RefSeq.col=TRUE)
        ans <- .format_Ensembl_chrom_info(seq_regions)
        .ENSEMBL_cached_chrom_info[[caching_key]] <- ans
    }
    if (assembled.molecules.only) {
        keep_idx <- which(ans[ , "coord_system.name"] %in% "chromosome")
        ans <- S4Vectors:::extract_data_frame_rows(ans, keep_idx)
    }
    if (!include.no_ref.sequences) {
        keep_idx <- which(!ans[ , "non_ref"])
        ans <- S4Vectors:::extract_data_frame_rows(ans, keep_idx)
    }
    ans
}

getChromInfoFromEnsembl <- function(dataset,
    release=NA, use.grch37=FALSE, kingdom=NA,
    assembled.molecules.only=FALSE,
    coord.system.names=c("chromosome", "scaffold"),
    include.no_ref.sequences=FALSE,
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
    if (!isTRUEorFALSE(include.no_ref.sequences))
        stop(wmsg("'include.no_ref.sequences' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(as.Seqinfo))
        stop(wmsg("'as.Seqinfo' must be TRUE or FALSE"))

    ans <- .get_chrom_info_for_unregistered_Ensembl_dataset(core_url,
                assembled.molecules.only=assembled.molecules.only,
                coord.system.names=coord.system.names,
                include.no_ref.sequences=include.no_ref.sequences,
                recache=recache)

    if (!as.Seqinfo)
        return(ans)
    Seqinfo(seqnames=ans[ , "name"],
            seqlengths=ans[ , "length"],
            isCircular=ans[ , "circular"],
            genome=genome)
}

