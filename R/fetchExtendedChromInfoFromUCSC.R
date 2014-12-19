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

### 'NCBI_seqlevel' and 'NCBI_accns' must be parallel vectors.
.map_UCSC_seqlevel_to_NCBI_seqlevel <- function(UCSC_seqlevel,
                                                NCBI_seqlevel,
                                                NCBI_accns,
                                                special_mappings=NULL)
{
    ans <- rep.int(NA_integer_, length(UCSC_seqlevel))

    ## 1. Handle special mappings.
    if (!is.null(special_mappings)) {
        m1 <- match(names(special_mappings), UCSC_seqlevel)
        if (any(is.na(m1)))
            stop("'special_mappings' has names not in 'UCSC_seqlevel'")
        m2 <- match(special_mappings, NCBI_seqlevel)
        if (any(is.na(m2)))
            stop("'special_mappings' has values not in 'NCBI_seqlevel'")
        ans[m1] <- m2
    }
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## From now on, matching is case insensitive.
    UCSC_seqlevel <- tolower(UCSC_seqlevel)
    NCBI_seqlevel <- tolower(NCBI_seqlevel)

    ## 2. We assign based on exact matching of the seqlevels.
    m <- match(UCSC_seqlevel[unmapped_idx], NCBI_seqlevel)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 3. We assign based on exact matching of the seqlevels (after
    ##    removal of the "chr" prefix).
    nochr_prefix_seqlevels <- sub("^chr", "", UCSC_seqlevel[unmapped_idx])
    m <- match(nochr_prefix_seqlevels, NCBI_seqlevel)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 4. We assign based on the number embedded in the UCSC chromosome name.
    split_seqlevels <- strsplit(UCSC_seqlevel[unmapped_idx], "_")
    nparts <- elementLengths(split_seqlevels)
    idx2 <- which(nparts >= 2L)
    if (length(idx2) != 0L) {
        ans[unmapped_idx[idx2]] <- sapply(idx2,
            function(i2) {
                part2 <- split_seqlevels[[i2]][2L]
                part2 <- sub("v", ".", part2, fixed=TRUE)
                m <- grep(part2, NCBI_accns, ignore.case=TRUE)
                if (length(m) >= 2L)
                    stop("cannot map ", UCSC_seqlevel[unmapped_idx[i2]],
                         " to a unique NCBI seqlevel")
                if (length(m) == 0L)
                    return(NA_integer_)
                m
            })
    }
    ans
}

standard_fetch_extended_ChromInfo_from_UCSC <- function(
                        genome,
                        circ_seqs,
                        refseq_assembly_accession,
                        AssemblyUnits,
                        special_mappings,
                        unmapped_seqs,
                        goldenPath_url,
                        quiet)
{
    chrominfo <- fetch_ChromInfo_from_UCSC(genome,
                                goldenPath_url=goldenPath_url)
    UCSC_seqlevel <- chrominfo[ , "chrom"]
    ans <- data.frame(UCSC_seqlevel=UCSC_seqlevel,
                      UCSC_seqlength=chrominfo[ , "size"],
                      circular=UCSC_seqlevel %in% circ_seqs,
                      stringsAsFactors=FALSE)
    if (is.null(refseq_assembly_accession)) {
        if (!quiet)
            warning(genome, " UCSC genome is not based on an NCBI assembly")
        oo <- order(rankSeqlevels(ans[ , "UCSC_seqlevel"]))
        ans <- ans[oo, , drop=FALSE]
        rownames(ans) <- NULL
        return(ans)
    }
    if (length(unmapped_seqs) != 0L) {
        unmapped_seqs_role <- rep.int(names(unmapped_seqs),
                                      elementLengths(unmapped_seqs))
        unmapped_seqs <- unlist(unmapped_seqs, use.names=FALSE)
        unmapped_idx <- match(unmapped_seqs, UCSC_seqlevel)
        stopifnot(!any(is.na(unmapped_idx)))
        if (!quiet)
            warning("NCBI seqlevel was set to NA for ", genome,
                    " UCSC seqlevel(s) not in the NCBI assembly: ",
                    paste(unmapped_seqs, collapse=", "))
    }
    assembly_report <- fetch_assembly_report(refseq_assembly_accession,
                                             AssemblyUnits=AssemblyUnits)
    #stopifnot(nrow(chrominfo) == nrow(assembly_report) + length(unmapped_seqs))
    NCBI_seqlevel <- assembly_report[ , "SequenceName"]
    GenBankAccn <- assembly_report[ , "GenBankAccn"]
    m <- .map_UCSC_seqlevel_to_NCBI_seqlevel(UCSC_seqlevel,
                            NCBI_seqlevel,
                            GenBankAccn,
                            special_mappings=special_mappings)
    if (length(unmapped_seqs) != 0L)
        stopifnot(all(is.na(m)[unmapped_idx]))
    unexpectedly_unmapped_idx <-
        which(is.na(m) & !(UCSC_seqlevel %in% unmapped_seqs))
    if (length(unexpectedly_unmapped_idx) != 0L)
        stop("cannot map ", genome, " UCSC seqlevel(s) ",
             paste(UCSC_seqlevel[unexpectedly_unmapped_idx], collapse=", "),
             " to an NCBI seqlevel")
    GenBankAccn[which(GenBankAccn == "na")] <- NA_character_
    SequenceRole <- factor(assembly_report[ , "SequenceRole"],
                           levels=c("assembled-molecule",
                                    "alt-scaffold",
                                    "unlocalized-scaffold",
                                    "unplaced-scaffold",
                                    "pseudo-scaffold"))
    stopifnot(identical(is.na(SequenceRole),
                        is.na(assembly_report[ , "SequenceRole"])))
    ans <- cbind(ans, NCBI_seqlevel=NCBI_seqlevel[m],
                      SequenceRole=SequenceRole[m],
                      GenBankAccn=GenBankAccn[m],
                      stringsAsFactors=FALSE)
    if (length(unmapped_seqs) != 0L)
        ans[unmapped_idx, "SequenceRole"] <- unmapped_seqs_role
    oo <- order(as.integer(ans[ , "SequenceRole"]),
                rankSeqlevels(ans[ , "UCSC_seqlevel"]))
    ans <- ans[oo, , drop=FALSE]
    rownames(ans) <- NULL
    ans
}

### See http://genome.ucsc.edu/FAQ/FAQreleases.html#release1 for the list
### of all UCSC genomes.
SUPPORTED_UCSC_GENOMES <- list(

### Human
    hg38=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000001405.26",
        special_mappings=c(chrM="MT")
    ),

    hg19=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000001405.13",
        ## Special renaming of the 9 alternate scaffolds:
        special_mappings=c(chr4_ctg9_hap1="HSCHR4_1_CTG9",
                           chr6_apd_hap1="HSCHR6_MHC_APD_CTG1",
                           chr6_cox_hap2="HSCHR6_MHC_COX_CTG1",
                           chr6_dbb_hap3="HSCHR6_MHC_DBB_CTG1",
                           chr6_mann_hap4="HSCHR6_MHC_MANN_CTG1",
                           chr6_mcf_hap5="HSCHR6_MHC_MCF_CTG1",
                           chr6_qbl_hap6="HSCHR6_MHC_QBL_CTG1",
                           chr6_ssto_hap7="HSCHR6_MHC_SSTO_CTG1",
                           chr17_ctg5_hap1="HSCHR17_1_CTG5"),
        unmapped_seqs=list(`assembled-molecule`="chrM")
    ),

    hg18=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000001405.12",
        special_mappings=c(chr6_cox_hap1="Hs6_111610_36",
                           chr22_h2_hap1="Hs22_111678_36"),
        unmapped_seqs=list(
            `assembled-molecule`="chrM",
            `pseudo-scaffold`=paste0("chr",
                c("5_h2_hap1", "6_qbl_hap2",
                  paste0(c((1:22)[-c(12, 14, 20)], "X"), "_random"))))
        ),

### Cow
    bosTau8=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000003055.5",
        special_mappings=c(chrM="MT")
    ),

    bosTau7=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000003205.5",
        unmapped_seqs=list(`assembled-molecule`="chrM")
    ),

    bosTau6=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000003055.4",
        special_mappings=c(chrM="MT")
    ),

### Dog
    canFam3=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000002285.3",
        special_mappings=c(chr1="chr01", chr2="chr02", chr3="chr03",
                           chr4="chr04", chr5="chr05", chr6="chr06",
                           chr7="chr07", chr8="chr08", chr9="chr09",
                           chrM="MT")
    ),

    canFam2=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM"
    ),

    canFam1=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM"
    ),

### Mouse
    mm10=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000001635.20",
        AssemblyUnits=c("C57BL/6J", "non-nuclear"),
        special_mappings=c(chrM="MT")
    ),

    mm9=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000001635.18",
        AssemblyUnits=c("C57BL/6J", "non-nuclear"),
        special_mappings=c(chrM="MT"),
        unmapped_seqs=list(
            `pseudo-scaffold`=paste0("chr",
                c(1, 3:5, 7:9, 13, 16:17, "X", "Y", "Un"), "_random"))
    ),

### Rat
    rn6=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000001895.5",
        special_mappings=c(chrM="MT")
    ),

### D. melanogaster
    dm6=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000001215.4",
        special_mappings=c(chrM="MT")
    ),

    dm3=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000001215.2",
        special_mappings=c(chrM="MT", chrU="Un"),
        unmapped_seqs=list(`pseudo-scaffold`="chrUextra")
    ),

### C. elegans
    ce10=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM"
    ),

    ce6=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000002985.1",
        special_mappings=c(chrM="MT")
    ),

    ce4=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM"
    ),

    ce2=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM"
    ),

### Yeast
    sacCer3=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        refseq_assembly_accession="GCF_000146045.2",
        special_mappings=c(chrM="MT")
    ),

    sacCer2=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs=c("chrM", "2micron")
    )
    ## more to come...
)

fetchExtendedChromInfoFromUCSC <- function(genome,
        goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
        quiet=FALSE)
{
    if (!isSingleString(genome))
        stop("'genome' must be a single string")
    if (!isTRUEorFALSE(quiet))
        stop("'quiet' must be TRUE or FALSE")
    idx <- match(genome, names(SUPPORTED_UCSC_GENOMES))
    if (is.na(idx))
        stop("genome \"", genome, "\" is not supported")
    supported_genome <- SUPPORTED_UCSC_GENOMES[[idx]]
    FUN <- get(supported_genome$FUN)
    FUN(genome=names(SUPPORTED_UCSC_GENOMES)[idx],
        circ_seqs=supported_genome$circ_seqs,
        refseq_assembly_accession=supported_genome$refseq_assembly_accession,
        AssemblyUnits=supported_genome$AssemblyUnits,
        special_mappings=supported_genome$special_mappings,
        unmapped_seqs=supported_genome$unmapped_seqs,
        goldenPath_url=goldenPath_url,
        quiet=quiet)
}

