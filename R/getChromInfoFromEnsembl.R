### =========================================================================
### getChromInfoFromEnsembl()
### -------------------------------------------------------------------------
###



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getChromInfoFromEnsembl()
###

.normarg_coord.system.names <- function(coord.system.names)
{
    if (is.null(coord.system.names))
        return(NULL)
    if (!is.character(coord.system.names))
        stop(wmsg("'coord.system.names' must be a character vector or NULL"))
    if (anyNA(coord.system.names) ||
        !all(nzchar(coord.system.names)) ||
        anyDuplicated(coord.system.names))
       stop(wmsg("'coord.system.names' cannot contain NAs, ",
                 "empty strings, or duplicates"))
    ## We sort to push normalization even further. This important because
    ## we will encode the normalized 'coord.system.names' vector in the
    ## caching key.
    sort(coord.system.names)
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
        ans <- fetch_seq_regions_from_Ensembl_ftp(core_url,
                                      coord_system_names=coord.system.names,
                                      add.toplevel.col=TRUE,
                                      add.non_ref.col=TRUE,
                                      add.INSDC.col=TRUE,
                                      add.RefSeq.col=TRUE)
        ## Add column "circular".
        circular <- make_circ_flags_from_circ_seqs(ans[ , "name"])
        ans <- cbind(ans, circular=circular)

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

