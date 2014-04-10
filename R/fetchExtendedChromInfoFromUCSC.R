### =========================================================================
### fetchExtendedChromInfoFromUCSC()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some low-level helpers.
###

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
         paste(UCSC_seqlevels[not_ok_idx], collapse=", "))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetchExtendedChromInfoFromUCSC()
###

.standard_fetch_extended_ChromInfo_from_UCSC <- function(genome,
                                                         refseq_assembly_id,
                                                         AssemblyUnits,
                                                         goldenPath_url)
{
    chrominfo <- fetch_ChromInfo_from_UCSC(genome,
                                goldenPath_url=goldenPath_url)
    UCSC_seqlevels <- chrominfo$chrom
    NCBI_accn2seqlevel <- fetch_GenBankAccn2seqlevel_from_NCBI(
                                refseq_assembly_id,
                                AssemblyUnits=AssemblyUnits)
    stopifnot(length(NCBI_accn2seqlevel) == length(UCSC_seqlevels))
    UCSC_accns <- .assign_GenBankAccns_to_UCSCseqlevels(
                                UCSC_seqlevels,
                                NCBI_accn2seqlevel)
    stopifnot(anyDuplicated(UCSC_accns) == 0L)
    NCBI_seqlevels <- unname(NCBI_accn2seqlevel[UCSC_accns])
    data.frame(UCSC_seqlevels=UCSC_seqlevels,
               NCBI_seqlevels=NCBI_seqlevels,
               accn=UCSC_accns,
               seqlengths=chrominfo$size,
               stringsAsFactors=FALSE)
}

.SUPPORTED_GENOMES <- list(

    hg38=list(FUN=.standard_fetch_extended_ChromInfo_from_UCSC,
              refseq_assembly_id="GCF_000001405.26"),

    mm10=list(FUN=.standard_fetch_extended_ChromInfo_from_UCSC,
              refseq_assembly_id="GCF_000001635.20",
              AssemblyUnits=c("C57BL/6J", "non-nuclear"))

    ## It's impossible to map the sequences in mm9 with the sequences in
    ## MGSCv37 because the mapping doesn't seem to be one-to-one. For example,
    ## it seems that the 102 unlocalized sequences on chromosome Y have been
    ## merged into a single sequence in mm9, the chrY_random sequence.
    ## So there is no way we can support mm9 :-/
    #mm9=list(FUN=.standard_fetch_extended_ChromInfo_from_UCSC,
    #          refseq_assembly_id="GCF_000001635.18",
    #          AssemblyUnits=c("C57BL/6J", "non-nuclear"))

    ## more to come...
)

.genome2idx <- function(genome)
{
    refseq_assembly_id <- lookup_refseq_assembly_id(genome)
    if (is.na(refseq_assembly_id))
        return(NA_integer_)
    refseq_assembly_ids <- sapply(.SUPPORTED_GENOMES,
                                  `[[`, "refseq_assembly_id",
                                  USE.NAMES=FALSE)
    match(refseq_assembly_id, refseq_assembly_ids)
}

### Only supports UCSC genomes. However 'genome' can be:
###   (a) a UCSC assembly name (e.g. "hg38");
###   (b) a RefSeq Assembly ID (e.g. "GCF_000001405.26");
###   (c) a GenBank Assembly ID (e.g. "GCA_000001405.15");
###   (d) an NCBI assembly name (e.g. "GRCh38").
fetchExtendedChromInfoFromUCSC <- function(genome,
        goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
{
    if (!.isSingleString(genome))
        stop("'genome' must be a single string")
    idx <- match(genome, names(.SUPPORTED_GENOMES))
    if (is.na(idx)) {
        idx <- .genome2idx(genome)
        if (is.na(idx))
            stop("genome \"", genome, "\" is not supported")
    }
    supported_genome <- .SUPPORTED_GENOMES[[idx]]
    supported_genome$FUN(names(.SUPPORTED_GENOMES)[idx],
                         supported_genome$refseq_assembly_id,
                         supported_genome$AssemblyUnits,
                         goldenPath_url)
}

