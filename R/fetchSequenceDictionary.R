### =========================================================================
### fetchSequenceDictionary()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .fetch_seq_dict_for_hg38()
###

fetch_GenBankAccn2seqlevel_for_GRCh38 <- function()
{
    NCBI_assembly_report <- fetch_assembly_report("GCF_000001405.26")
    GenBank_accn <- NCBI_assembly_report[[5L]]
    stopifnot(!any(duplicated(GenBank_accn)))
    ans <- NCBI_assembly_report[[1L]]
    names(ans) <- GenBank_accn
    ans
}

.fetch_chrominfo_for_hg38 <- function()
{
    hg38_chrominfo_url <- paste(
        "http://hgdownload.soe.ucsc.edu/goldenPath",
        "hg38/database/chromInfo.txt.gz", sep="/")
    destfile <- tempfile()
    download.file(hg38_chrominfo_url, destfile, quiet=TRUE)
    colnames <- c("chrom", "size", "fileName")
    read.table(destfile, sep="\t", quote="",
                         col.names=colnames, comment.char="",
                         stringsAsFactors=FALSE)
}

.make_GenBankAccn2seqlevel_for_hg38 <- function(hg38_chrominfo,
                                                GRCh38_accn2seqlevel)
{
    suppressMessages(require(IRanges, quietly=TRUE))  # for elementLengths()
    hg38_seqlevels <- hg38_chrominfo$chrom
    tmp <- lapply(names(GRCh38_accn2seqlevel), grep, hg38_seqlevels)
    idx1 <- which(elementLengths(tmp) == 1L)
    tmp[idx1] <- as.list(hg38_seqlevels[unlist(tmp[idx1], use.names=FALSE)])
    NCBI_main_chromosomes <- c(1:22, "X", "Y", "MT")
    UCSC_main_chromosomes <- paste0("chr", c(1:22, "X", "Y", "M"))
    idx0 <- which(elementLengths(tmp) == 0L)
    stopifnot(all(GRCh38_accn2seqlevel[idx0] == NCBI_main_chromosomes))
    tmp[idx0] <- as.list(UCSC_main_chromosomes)
    stopifnot(all(elementLengths(tmp) == 1L))
    ans <- unlist(tmp, use.names=FALSE)
    names(ans) <- names(GRCh38_accn2seqlevel)
    ans
}

.fetch_seq_dict_for_hg38 <- function()
{
    GRCh38_accn2seqlevel <- fetch_GenBankAccn2seqlevel_for_GRCh38()
    hg38_chrominfo <- .fetch_chrominfo_for_hg38()
    hg38_accn2seqlevel <- .make_GenBankAccn2seqlevel_for_hg38(
                              hg38_chrominfo,
                              GRCh38_accn2seqlevel)
    GRCh38_seqlevels <- unname(GRCh38_accn2seqlevel)
    hg38_seqlevels <- unname(hg38_accn2seqlevel)
    seqlengths <- hg38_chrominfo$size[match(hg38_seqlevels,
                                            hg38_chrominfo$chrom)]
    data.frame(accn=names(GRCh38_accn2seqlevel),
               GRCh38_seqlevels=GRCh38_seqlevels,
               hg38_seqlevels=hg38_seqlevels,
               seqlengths=seqlengths,
               stringsAsFactors=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetchSequenceDictionary()
###

.SEQUENCE_DICTIONARIES <- list(
    hg38=list(.fetch_seq_dict_for_hg38, "GCF_000001405.26")
    ## more to come...
)

.assembly2idx <- function(assembly)
{
    refseq_assembly_id <- lookup_refseq_assembly_id(assembly)
    if (is.na(refseq_assembly_id))
        return(NA_integer_)
    refseq_assembly_ids <- sapply(.SEQUENCE_DICTIONARIES, `[[`, 2L,
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
    FUN <- .SEQUENCE_DICTIONARIES[[idx]][[1L]]
    FUN()
}

