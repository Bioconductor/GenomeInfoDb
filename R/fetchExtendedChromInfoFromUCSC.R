### =========================================================================
### fetchExtendedChromInfoFromUCSC()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some low-level helpers.
###

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
    assembly_report <- fetch_assembly_report(assembly,
                                             AssemblyUnits=AssemblyUnits)
    GenBank_accn <- assembly_report[["GenBankAccn"]]
    if ("na" %in% GenBank_accn)
        stop("GenBankAccn field in assembly report for ",
             "\"", assembly, "\" contains \"na\"")
    stopifnot(anyDuplicated(GenBank_accn) == 0L)
    ans <- assembly_report[["SequenceName"]]
    names(ans) <- GenBank_accn
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetchExtendedChromInfoFromUCSC()
###

### 'NCBI_seqlevels' and 'NCBI_accns' must be parallel vectors.
.map_UCSC_seqlevels_to_NCBI_seqlevels <- function(UCSC_seqlevels,
                                                  NCBI_seqlevels,
                                                  NCBI_accns,
                                                  special_renamings=NULL)
{
    ans <- rep.int(NA_integer_, length(UCSC_seqlevels))

    ## 1. Handle special renamings.
    if (!is.null(special_renamings)) {
        m1 <- match(names(special_renamings), UCSC_seqlevels)
        if (any(is.na(m1)))
            stop("'special_renamings' has names not in 'UCSC_seqlevels'")
        m2 <- match(special_renamings, NCBI_seqlevels)
        if (any(is.na(m2)))
            stop("'special_renamings' has values not in 'NCBI_seqlevels'")
        ans[m1] <- m2
    }
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 2. We assign based on exact matching of the sequence names (after
    ##    removal of the "chr" prefix).
    nochr_prefix_seqlevels <- sub("^chr", "", UCSC_seqlevels[unmapped_idx])
    m <- match(nochr_prefix_seqlevels, NCBI_seqlevels)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 3. We assign based on the number embedded in the UCSC chromosome name.
    split_seqlevels <- strsplit(UCSC_seqlevels[unmapped_idx], "_")
    nparts <- elementLengths(split_seqlevels)
    idx2 <- which(nparts >= 2L)
    if (length(idx2) != 0L) {
        ans[unmapped_idx[idx2]] <- sapply(idx2,
            function(i2) {
                part2 <- split_seqlevels[[i2]][2L]
                part2 <- sub("v", ".", part2, fixed=TRUE)
                m <- grep(part2, NCBI_accns, ignore.case=TRUE)
                if (length(m) >= 2L)
                    stop("cannot map ", UCSC_seqlevels[unmapped_idx[i2]],
                         " to a unique NCBI seqlevel")
                if (length(m) == 0L)
                    return(NA_integer_)
                m
            })
    }
    ans
}

standard_fetch_extended_ChromInfo_from_UCSC <- function(genome,
                                                        circular,
                                                        refseq_assembly_id,
                                                        AssemblyUnits,
                                                        special_renamings,
                                                        unmapped,
                                                        goldenPath_url)
{
    chrominfo <- fetch_ChromInfo_from_UCSC(genome,
                                goldenPath_url=goldenPath_url)
    UCSC_seqlevels <- chrominfo$chrom
    ans <- data.frame(UCSC_seqlevels=UCSC_seqlevels,
                      UCSC_seqlengths=chrominfo$size,
                      circular=UCSC_seqlevels %in% circular,
                      stringsAsFactors=FALSE)
    if (is.null(refseq_assembly_id)) {
        warning(genome, " UCSC genome is not based on an NCBI assembly")
        return(ans)
    }
    if (!is.null(unmapped)) {
        unmapped_idx <- match(unmapped, UCSC_seqlevels)
        stopifnot(!any(is.na(unmapped_idx)))
        warning("NCBI seqlevel was set to NA for ", genome,
                " UCSC seqlevel(s) not in the NCBI assembly: ",
                paste(unmapped, collapse=", "))
    }
    assembly_report <- fetch_assembly_report(refseq_assembly_id,
                                             AssemblyUnits=AssemblyUnits)
    #stopifnot(nrow(chrominfo) == nrow(assembly_report) + length(unmapped))
    NCBI_seqlevels <- assembly_report$SequenceName
    GenBankAccn <- assembly_report$GenBankAccn
    m <- .map_UCSC_seqlevels_to_NCBI_seqlevels(UCSC_seqlevels,
                                NCBI_seqlevels,
                                GenBankAccn,
                                special_renamings=special_renamings)
    if (!is.null(unmapped))
        stopifnot(all(is.na(m)[unmapped_idx]))
    unexpectedly_unmapped_idx <-
        which(is.na(m) & !(UCSC_seqlevels %in% unmapped))
    if (length(unexpectedly_unmapped_idx) != 0L)
        stop("cannot map ", genome, " UCSC seqlevel(s) ",
             paste(UCSC_seqlevels[unexpectedly_unmapped_idx], collapse=", "),
             " to an NCBI seqlevel")
    GenBankAccn[which(GenBankAccn == "na")] <- NA_character_
    SequenceRole <- factor(assembly_report$SequenceRole,
                           levels=c("assembled-molecule",
                                    "unlocalized-scaffold",
                                    "unplaced-scaffold",
                                    "alt-scaffold",
                                    "pseudo-scaffold"))
    stopifnot(identical(is.na(SequenceRole),
                        is.na(assembly_report$SequenceRole)))
    ans <- cbind(ans, NCBI_seqlevels=NCBI_seqlevels[m],
                      SequenceRole=SequenceRole[m],
                      GenBankAccn=GenBankAccn[m],
                      stringsAsFactors=FALSE)
    ans[order(ans$SequenceRole), , drop=FALSE]
}

SUPPORTED_UCSC_GENOMES <- list(
    hg38=
        list(FUN="standard_fetch_extended_ChromInfo_from_UCSC",
             circular="chrM",
             refseq_assembly_id="GCF_000001405.26",
             special_renamings=c(chrM="MT")),
    hg19=
        list(FUN="standard_fetch_extended_ChromInfo_from_UCSC",
             circular="chrM",
             refseq_assembly_id="GCF_000001405.13",
             ## Special renaming of the 9 alternate scaffolds:
             special_renamings=c(chr6_apd_hap1="HSCHR6_MHC_APD_CTG1",
                                 chr6_cox_hap2="HSCHR6_MHC_COX_CTG1",
                                 chr6_dbb_hap3="HSCHR6_MHC_DBB_CTG1",
                                 chr6_mann_hap4="HSCHR6_MHC_MANN_CTG1",
                                 chr6_mcf_hap5="HSCHR6_MHC_MCF_CTG1",
                                 chr6_qbl_hap6="HSCHR6_MHC_QBL_CTG1",
                                 chr6_ssto_hap7="HSCHR6_MHC_SSTO_CTG1",
                                 chr4_ctg9_hap1="HSCHR4_1_CTG9",
                                 chr17_ctg5_hap1="HSCHR17_1_CTG5"),
             unmapped="chrM"),
    hg18=
        list(FUN="standard_fetch_extended_ChromInfo_from_UCSC",
             circular="chrM",
             refseq_assembly_id="GCF_000001405.12",
             AssemblyUnits="Primary Assembly",
             unmapped=paste0("chr",
                 c("M", paste0(c((1:22)[-c(12, 14, 20)], "X"), "_random"),
                   "5_h2_hap1", "6_cox_hap1", "6_qbl_hap2", "22_h2_hap1"))),
    mm10=
        list(FUN="standard_fetch_extended_ChromInfo_from_UCSC",
             circular="chrM",
             refseq_assembly_id="GCF_000001635.20",
             AssemblyUnits=c("C57BL/6J", "non-nuclear"),
             special_renamings=c(chrM="MT")),
    mm9=
        list(FUN="standard_fetch_extended_ChromInfo_from_UCSC",
             circular="chrM",
             refseq_assembly_id="GCF_000001635.18",
             AssemblyUnits=c("C57BL/6J", "non-nuclear"),
             special_renamings=c(chrM="MT"),
             unmapped=paste0("chr",
                 paste0(c(1, 3:5, 7:9, 13, 16:17, "X", "Y", "Un"), "_random"))),
    rn6=
        list(FUN="standard_fetch_extended_ChromInfo_from_UCSC",
             circular="chrM",
             refseq_assembly_id="GCA_000001895.4",
             special_renamings=c(chrM="MT")),
    dm6=
        list(FUN="standard_fetch_extended_ChromInfo_from_UCSC",
             circular="chrM",
             refseq_assembly_id="GCF_000001215.4",
             special_renamings=c(chrM="MT")),
    dm3=
        list(FUN="standard_fetch_extended_ChromInfo_from_UCSC",
             circular="chrM",
             refseq_assembly_id="GCF_000001215.2",
             special_renamings=c(chrM="MT", chrU="Un"),
             unmapped="chrUextra"),
    sacCer3=
        list(FUN="standard_fetch_extended_ChromInfo_from_UCSC",
             circular="chrM",
             refseq_assembly_id="GCF_000146045.2",
             special_renamings=c(chrM="MT")),
    sacCer2=
        list(FUN="standard_fetch_extended_ChromInfo_from_UCSC",
             circular=c("chrM", "2micron"))
    ## more to come...
)

.genome2idx <- function(genome)
{
    refseq_assembly_id <- lookup_refseq_assembly_id(genome)
    if (is.na(refseq_assembly_id))
        return(NA_integer_)
    refseq_assembly_ids <- sapply(SUPPORTED_UCSC_GENOMES,
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
    if (!isSingleString(genome))
        stop("'genome' must be a single string")
    idx <- match(genome, names(SUPPORTED_UCSC_GENOMES))
    if (is.na(idx)) {
        idx <- .genome2idx(genome)
        if (is.na(idx))
            stop("genome \"", genome, "\" is not supported")
    }
    supported_genome <- SUPPORTED_UCSC_GENOMES[[idx]]
    FUN <- get(supported_genome$FUN)
    FUN(genome=names(SUPPORTED_UCSC_GENOMES)[idx],
        circular=supported_genome$circular,
        refseq_assembly_id=supported_genome$refseq_assembly_id,
        AssemblyUnits=supported_genome$AssemblyUnits,
        special_renamings=supported_genome$special_renamings,
        unmapped=supported_genome$unmapped,
        goldenPath_url=goldenPath_url)
}

