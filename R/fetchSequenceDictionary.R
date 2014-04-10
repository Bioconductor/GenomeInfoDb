### =========================================================================
### fetchSequenceDictionary()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some low-level helpers.
###

### Used in BSgenome!
fetch_GenBankAccn2seqlevel_from_NCBI <- function(assembly, AssemblyUnits=NULL)
{
    assembly_report <- fetch_assembly_report(assembly)
    if (!is.null(AssemblyUnits)) {
        stopifnot(all(AssemblyUnits %in% assembly_report$AssemblyUnit))
        idx <- which(assembly_report$AssemblyUnit %in% AssemblyUnits)
        assembly_report <- assembly_report[idx, , drop=FALSE]
    }
    GenBank_accn <- assembly_report[["GenBankAccn"]]
    stopifnot(anyDuplicated(GenBank_accn) == 0L)
    ans <- assembly_report[["SequenceName"]]
    names(ans) <- GenBank_accn
    ans
}

### Use this in GenomicFeatures/R/makeTranscriptDbFromUCSC.R instead of
### internal utility .downloadChromInfoFromUCSC().
fetch_ChromInfo_from_UCSC <- function(genome,
        goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
{
    url <- paste(goldenPath_url, genome, "database/chromInfo.txt.gz", sep="/")
    destfile <- tempfile()
    download.file(url, destfile, quiet=TRUE)
    colnames <- c("chrom", "size", "fileName")
    read.table(destfile, sep="\t", quote="",
                         col.names=colnames, comment.char="",
                         stringsAsFactors=FALSE)
}

### WARNING! Using this for hg18 will assign the *wrong* GenBank accession
### number to chrM!
.assign_GenBankAccns_to_UCSCseqlevels <- function(UCSC_seqlevels,
                                                  NCBI_accn2seqlevel)
{
    suppressMessages(require(IRanges, quietly=TRUE))  # for elementLengths()
    ans <- rep.int(NA_character_, length(UCSC_seqlevels))

    ## 1. We assign based on exact matching of the sequence names (after
    ##    removal of the "chr" prefix).
    nochr_prefix_seqlevels <- sub("^chr", "", UCSC_seqlevels)
    m <- match(nochr_prefix_seqlevels, NCBI_accn2seqlevel)
    ok_idx <- which(!is.na(m))
    ans[ok_idx] <- names(NCBI_accn2seqlevel)[m[ok_idx]]
    if (length(ok_idx) == length(ans))
        return(ans)

    ## 2. We assign based on the number embedded in the UCSC chromosome name.
    not_ok_idx <- which(is.na(ans))
    split_seqlevels <- strsplit(UCSC_seqlevels[not_ok_idx], "_")
    nparts <- elementLengths(split_seqlevels)
    idx2 <- which(nparts >= 2L)
    ans[not_ok_idx[idx2]] <- sapply(idx2,
        function(i2) {
            part2 <- split_seqlevels[[i2]][2L]
            part2 <- sub("v", ".", part2, fixed=TRUE)
            accn <- grep(part2, names(NCBI_accn2seqlevel),
                         value=TRUE, fixed=TRUE)
            if (length(accn) >= 2L)
                stop("cannot assign a GenBank accession number to ",
                     UCSC_seqlevels[not_ok_idx[i2]])
            if (length(accn) == 0L)
                return(NA_character_)
            accn
        })
    not_ok_idx <- which(is.na(ans))
    if (length(not_ok_idx) == 0L)
        return(ans)

    ## 3. Handle common sequence renamings.
    common_renamings <- c(
        chrM="MT"
        ## add more common renamings...
    )
    m <- match(names(common_renamings), UCSC_seqlevels)
    ok_idx <- which(!is.na(m))
    m2 <- match(common_renamings[ok_idx], NCBI_accn2seqlevel)
    ans[m[ok_idx]] <- names(NCBI_accn2seqlevel)[m2]
    not_ok_idx <- which(is.na(ans))
    if (length(not_ok_idx) == 0L)
        return(ans)

    ## 4. Fail.
    stop("failed to assign a GenBank accession number to: ",
         paste(UCSC_seqlevels[not_ok_idx], sep=", "))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetchSequenceDictionary()
###

.fetch_seq_dict_for_hg38 <- function()
{
    GRCh38_accn2seqlevel <- fetch_GenBankAccn2seqlevel_from_NCBI(
                                "GCF_000001405.26")
    hg38_chrominfo <- fetch_ChromInfo_from_UCSC("hg38")
    hg38_seqlevels <- hg38_chrominfo$chrom
    stopifnot(length(hg38_seqlevels) == length(GRCh38_accn2seqlevel))
    hg38_accns <- .assign_GenBankAccns_to_UCSCseqlevels(
                                hg38_seqlevels,
                                GRCh38_accn2seqlevel)
    stopifnot(anyDuplicated(hg38_accns) == 0L)

    GRCh38_seqlevels <- unname(GRCh38_accn2seqlevel[hg38_accns])
    seqlengths <- hg38_chrominfo$size[match(hg38_seqlevels,
                                            hg38_chrominfo$chrom)]
    data.frame(accn=hg38_accns,
               GRCh38_seqlevels=GRCh38_seqlevels,
               hg38_seqlevels=hg38_seqlevels,
               seqlengths=seqlengths,
               stringsAsFactors=FALSE)
}

.fetch_seq_dict_for_mm10 <- function()
{
    GRCm38_accn2seqlevel <- fetch_GenBankAccn2seqlevel_from_NCBI(
                                "GCF_000001635.20",
                                AssemblyUnits=c("C57BL/6J", "non-nuclear"))
    mm10_chrominfo <- fetch_ChromInfo_from_UCSC("mm10")
    mm10_seqlevels <- mm10_chrominfo$chrom
    stopifnot(length(mm10_seqlevels) == length(GRCm38_accn2seqlevel))
    mm10_accns <- .assign_GenBankAccns_to_UCSCseqlevels(
                                mm10_seqlevels,
                                GRCm38_accn2seqlevel)
    stopifnot(anyDuplicated(names(mm10_seqlevels)) == 0L)

    GRCm38_seqlevels <- unname(GRCm38_accn2seqlevel[mm10_accns])
    seqlengths <- mm10_chrominfo$size[match(mm10_seqlevels,
                                            mm10_chrominfo$chrom)]
    data.frame(accn=mm10_accns,
               GRCm38_seqlevels=GRCm38_seqlevels,
               mm10_seqlevels=mm10_seqlevels,
               seqlengths=seqlengths,
               stringsAsFactors=FALSE)
}

.SEQUENCE_DICTIONARIES <- list(
    hg38=list(FUN=.fetch_seq_dict_for_hg38,
              refseq_assembly_id="GCF_000001405.26"),
    mm10=list(FUN=.fetch_seq_dict_for_mm10,
              refseq_assembly_id="GCF_000001635.20",
              AssemblyUnits=c("C57BL/6J", "non-nuclear"))
    ## more to come...
)

.assembly2idx <- function(assembly)
{
    refseq_assembly_id <- lookup_refseq_assembly_id(assembly)
    if (is.na(refseq_assembly_id))
        return(NA_integer_)
    refseq_assembly_ids <- sapply(.SEQUENCE_DICTIONARIES,
                                  `[[`, "refseq_assembly_id",
                                  USE.NAMES=FALSE)
    match(refseq_assembly_id, refseq_assembly_ids)
}

### Only supports UCSC assemblies. However 'assembly' can be:
###   (a) a UCSC assembly name (e.g. "hg38");
###   (b) a RefSeq Assembly ID (e.g. "GCF_000001405.26");
###   (c) a GenBank Assembly ID (e.g. "GCA_000001405.15");
###   (d) an NCBI assembly name (e.g. "GRCh38").
fetchSequenceDictionary <- function(assembly)
{
    if (!.isSingleString(assembly))
        stop("'assembly' must be a single string")
    idx <- match(assembly, names(.SEQUENCE_DICTIONARIES))
    if (is.na(idx)) {
        idx <- .assembly2idx(assembly)
        if (is.na(idx))
            stop("assembly \"", assembly, "\" is not supported")
    }
    FUN <- .SEQUENCE_DICTIONARIES[[idx]][["FUN"]]
    FUN()
}

