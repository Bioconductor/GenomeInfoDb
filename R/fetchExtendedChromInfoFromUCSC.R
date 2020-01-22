### =========================================================================
### fetchExtendedChromInfoFromUCSC()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some low-level helpers.
###

.fetch_table_from_UCSC <- function(genome, table, colnames=NULL,
        goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
{
    url <- paste(goldenPath_url, genome,
                 "database", paste0(table, ".txt.gz"), sep="/")
    destfile <- tempfile()
    download.file(url, destfile, quiet=TRUE)
    if (is.null(colnames)) {
        read.table(destfile, sep="\t", quote="",
                             comment.char="",
                             stringsAsFactors=FALSE)
    } else {
        read.table(destfile, sep="\t", quote="",
                             col.names=colnames,
                             comment.char="",
                             stringsAsFactors=FALSE)
    }
}

### NOT exported but used in the GenomicFeatures package.
fetch_ChromInfo_from_UCSC <- function(genome,
    goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
{
    colnames <- c("chrom", "size", "fileName")
    .fetch_table_from_UCSC(genome, "chromInfo", colnames,
                           goldenPath_url=goldenPath_url)
}

### NOT exported but used in the BSgenome package.
fetch_GenBankAccn2seqlevel_from_NCBI <- function(assembly, AssemblyUnits=NULL)
{
    assembly_report <- fetch_assembly_report(assembly,
                                             AssemblyUnits=AssemblyUnits)
    GenBank_accn <- assembly_report[["GenBankAccn"]]
    if ("na" %in% GenBank_accn)
        stop(wmsg("GenBankAccn field in assembly report for ",
                  "\"", assembly, "\" contains \"na\""))
    stopifnot(anyDuplicated(GenBank_accn) == 0L)
    ans <- assembly_report[["SequenceName"]]
    names(ans) <- GenBank_accn
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetchExtendedChromInfoFromUCSC()
###

..safe_match <- function(query, NCBI_accn)
{
    hits <- findMatches(query, NCBI_accn)
    q_hits <- queryHits(hits)
    s_hits <- subjectHits(hits)
    ambig_q_idx <- which(duplicated(q_hits))
    if (length(ambig_q_idx) != 0L) {
        ambig_idx <- unique(q_hits[ambig_q_idx])
        in1string <- paste0(names(query)[ambig_idx], collapse=", ")
        stop(wmsg("UCSC seqlevel(s) matching more than one accession number: ",
                  in1string))
    }
    ambig_s_idx <- which(duplicated(s_hits))
    if (length(ambig_s_idx) != 0L) {
        ambig_idx <- unique(s_hits[ambig_s_idx])
        in1string <- paste0(NCBI_accn[ambig_idx], collapse=", ")
        stop(wmsg("accession number(s) with more than one UCSC seqlevel ",
                  "match: ", in1string))
    }
    ans <- rep.int(NA_integer_, length(query))
    ans[q_hits] <- s_hits
    ans
}

.match_UCSC_seqlevel_to_NCBI_accn <- function(UCSC_seqlevel, NCBI_accn,
                                              accn_suffix="")
{
    query <- sub("-", ".", UCSC_seqlevel, fixed=TRUE)
    query <- paste0(query, accn_suffix)
    names(query) <- UCSC_seqlevel
    ..safe_match(tolower(query), tolower(NCBI_accn))
}

.match_UCSC_seqlevel_part2_to_NCBI_accn <- function(UCSC_seqlevel, NCBI_accn,
                                                    accn_prefix="")
{
    ans <- rep.int(NA_integer_, length(UCSC_seqlevel))
    seqlevel_parts <- strsplit(UCSC_seqlevel, "_")
    nparts <- elementNROWS(seqlevel_parts)
    idx2 <- which(nparts >= 2L)
    if (length(idx2) == 0L)
        return(ans)
    offsets <- c(0L, cumsum(nparts[idx2[-length(idx2)]]))
    query <- unlist(seqlevel_parts[idx2], use.names=FALSE)[offsets + 2L]
    query <- sub("v", ".", query, fixed=TRUE)
    unversioned_idx <- grep(".", query, fixed=TRUE, invert=TRUE)
    if (length(unversioned_idx) != 0L)
        query[unversioned_idx] <- paste0(query[unversioned_idx], ".1")
    query <- paste0(accn_prefix, query)
    names(query) <- UCSC_seqlevel[idx2]
    ans[idx2] <- ..safe_match(tolower(query), tolower(NCBI_accn))
    ans
}

### 'NCBI_seqlevel' and 'NCBI_accn' must be parallel vectors.
.map_UCSC_seqlevel_to_NCBI_seqlevel <- function(UCSC_seqlevel,
                                                NCBI_seqlevel,
                                                NCBI_accn,
                                                special_mappings=NULL)
{
    ans <- rep.int(NA_integer_, length(UCSC_seqlevel))

    ## 1. Handle special mappings.
    if (!is.null(special_mappings)) {
        m1 <- match(names(special_mappings), UCSC_seqlevel)
        if (any(is.na(m1)))
            stop(wmsg("'special_mappings' contains sequence names ",
                      "not in 'UCSC_seqlevel'"))
        m2 <- match(special_mappings, NCBI_seqlevel)
        if (any(is.na(m2)))
            stop(wmsg("'special_mappings' has values not in 'NCBI_seqlevel'"))
        ans[m1] <- m2
    }
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 2. We assign based on exact matching (case insensitive) of the
    ##    seqlevels.
    ucsc_seqlevel <- tolower(UCSC_seqlevel)
    ncbi_seqlevel <- tolower(NCBI_seqlevel)
    m <- match(ucsc_seqlevel[unmapped_idx], ncbi_seqlevel)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 3. We assign based on exact matching (case insensitive) of the
    ##    seqlevels after removal of the "chr" prefix.
    nochr_ucsc_seqlevel <- sub("^chr", "", ucsc_seqlevel[unmapped_idx])
    nochr_ncbi_seqlevel <- sub("^chr", "", ncbi_seqlevel)
    m <- match(nochr_ucsc_seqlevel, nochr_ncbi_seqlevel)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 4. We assign based on accession number found in UCSC seqlevel.
    m <- .match_UCSC_seqlevel_to_NCBI_accn(UCSC_seqlevel[unmapped_idx],
                                           NCBI_accn)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 5. We assign based on accession number found in UCSC seqlevel after
    ##    adding .1 suffix to it.
    m <- .match_UCSC_seqlevel_to_NCBI_accn(UCSC_seqlevel[unmapped_idx],
                                           NCBI_accn, accn_suffix=".1")
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 6. We assign based on accession number found in part 2 of UCSC seqlevel.
    m <- .match_UCSC_seqlevel_part2_to_NCBI_accn(UCSC_seqlevel[unmapped_idx],
                                                 NCBI_accn)
    ok_idx <- which(!is.na(m))
    ans[unmapped_idx[ok_idx]] <- m[ok_idx]
    unmapped_idx <- which(is.na(ans))
    if (length(unmapped_idx) == 0L)
        return(ans)

    ## 7. We assign based on accession number found in part 2 of UCSC seqlevel
    ##    after adding AAD prefix to it.
    ans[unmapped_idx] <- .match_UCSC_seqlevel_part2_to_NCBI_accn(
                                    UCSC_seqlevel[unmapped_idx], NCBI_accn,
                                    accn_prefix="AAD")
    ans
}

.fetch_Ensembl_seqlevels_from_UCSC <-
    function(genome, Ensembl_mapping_tables, UCSC_seqlevels, duplicated_contigs,
             goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
{
    stopifnot(identical(Ensembl_mapping_tables, "chromAlias") ||
              identical(Ensembl_mapping_tables,
                        c("chromAlias", "ucscToEnsembl")))
    ensembl_seqlevels <- rep(NA_character_, length(UCSC_seqlevels))
    names(ensembl_seqlevels) <- UCSC_seqlevels

    ## 1. Fill in what is already provided by the package GenomeInfoDb.
    org <- suppressWarnings(mapGenomeBuilds(genome)[1L, "organism"])
    if (!is.null(org)) {
        ## org2 <- .normalize_organism(org)
        ## available_orgs <- names(.getNamedFiles())
        ## if (org2 %in% available_orgs) {
        chroms <- try(genomeStyles(org), silent=TRUE)
        if (!is(chroms, "try-error")) {
            chroms <- chroms[chroms$UCSC %in% names(ensembl_seqlevels), ]
            ensembl_seqlevels[chroms$UCSC] <- chroms$Ensembl
        }
    }

    ## 2. Fetch info from UCSC.
    for (table in Ensembl_mapping_tables) {
        uscs_seqlevels_without_ensembl <-
            names(ensembl_seqlevels)[which(is.na(ensembl_seqlevels))]
        if (length(uscs_seqlevels_without_ensembl)) {
            if (table == "chromAlias") {
                colnames <- c("alias", "chrom", "source")
            } else {
                ## 'table' is "ucscToEnsembl"
                colnames <- c("chrom", "alias")
            }
            chroms <- .fetch_table_from_UCSC(genome, table, colnames,
                                             goldenPath_url=goldenPath_url)
            if ("source" %in% colnames(chroms)) {
                idx <- grep("ensembl", chroms[ , "source"])
                chroms <- chroms[idx, c("chrom", "alias")]
            }
            idx <- which(!(chroms[ , "alias"] %in% duplicated_contigs))
            chroms <- chroms[idx, ]
            rownames(chroms) <- chroms[ , "chrom"]
            keep_me <- uscs_seqlevels_without_ensembl %in% rownames(chroms)
            uscs_seqlevels_without_ensembl <-
                uscs_seqlevels_without_ensembl[keep_me]
            ensembl_seqlevels[uscs_seqlevels_without_ensembl] <-
                chroms[uscs_seqlevels_without_ensembl, "alias"]
        }
    }
    as.character(ensembl_seqlevels)
}

standard_fetch_extended_ChromInfo_from_UCSC <- function(
                                  genome,
                                  circ_seqs,
                                  assembly_accession,
                                  AssemblyUnits,
                                  special_mappings,
                                  unmapped_seqs,
                                  drop_unmapped,
                                  Ensembl_mapping_tables,
                                  duplicated_contigs,
                                  goldenPath_url,
                                  quiet)
{
    chrominfo <- fetch_ChromInfo_from_UCSC(genome,
                                           goldenPath_url=goldenPath_url)
    UCSC_seqlevel <- chrominfo[ , "chrom"]
    circular_idx <- match(circ_seqs, UCSC_seqlevel)
    if (any(is.na(circular_idx)))
        stop(wmsg("'circ_seqs' contains sequence names not in ",
                  genome, " genome"))
    circular <- logical(length(UCSC_seqlevel))
    circular[circular_idx] <- TRUE
    ans <- data.frame(UCSC_seqlevel=UCSC_seqlevel,
                      UCSC_seqlength=chrominfo[ , "size"],
                      circular=circular,
                      stringsAsFactors=FALSE)
    if (is.null(assembly_accession)) {
        if (!quiet)
            warning(wmsg(genome, " UCSC genome is not based on ",
                         "an NCBI assembly"))

        oo <- order(rankSeqlevels(ans[ , "UCSC_seqlevel"]))
        ans <- ans[oo, , drop=FALSE]
        rownames(ans) <- NULL
        return(ans)
    }
    if (length(unmapped_seqs) != 0L) {
        unmapped_seqs_role <- rep.int(names(unmapped_seqs),
                                      elementNROWS(unmapped_seqs))
        unmapped_seqs <- unlist(unmapped_seqs, use.names=FALSE)
        unmapped_idx <- match(unmapped_seqs, UCSC_seqlevel)
        stopifnot(!any(is.na(unmapped_idx)))
        if (!quiet)
            warning(wmsg("NCBI seqlevel was set to NA for ", genome,
                         " UCSC seqlevel(s) not in the NCBI assembly: ",
                         paste0(unmapped_seqs, collapse=", ")))
    }
    assembly_report <- fetch_assembly_report(assembly_accession,
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
    if (isTRUE(drop_unmapped)) {
        keep_idx <- which(!is.na(m))
        ans <- S4Vectors:::extract_data_frame_rows(ans, keep_idx)
        m <- m[keep_idx]
        ## Before we drop the unmapped UCSC seqlevels, we want to make sure
        ## that all the NCBI seqlevels are reverse mapped. For example, in
        ## the case of hg38, the chromInfo table at UCSC contains sequences
        ## that belong to GRCh38.p12 but not to GRCh38. So we want to make
        ## sure that after we drop these "foreign sequences", we are left
        ## with a one-to-one mapping between UCSC seqlevels and the 455 NCBI
        ## seqlevels that are in the assembly report for GRCh38.
        idx <- setdiff(seq_along(NCBI_seqlevel), m)
        if (length(idx) != 0L) {
            in1string <- paste0(NCBI_seqlevel[idx], collapse=", ")
            stop(wmsg("no UCSC seqlevel could be mapped to the following ",
                      "NCBI seqlevel(s): ", in1string))
        }
    } else {
        unexpectedly_unmapped_idx <-
            which(is.na(m) & !(UCSC_seqlevel %in% unmapped_seqs))
        if (length(unexpectedly_unmapped_idx) != 0L) {
            in1string <- paste0(UCSC_seqlevel[unexpectedly_unmapped_idx],
                                collapse=", ")
            stop(wmsg("cannot map the following UCSC seqlevel(s) to an ",
                      "NCBI seqlevel: ", in1string))
        }
    }
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
    if (length(unmapped_seqs) != 0L) {
        ## The order hardcoded in 'unmapped_seqs' is authority.
        unmapped_idx <- match(unmapped_seqs, ans[ , "UCSC_seqlevel"])
        ans[sort(unmapped_idx), ] <- ans[unmapped_idx, , drop=FALSE]
    }
    if (!is.null(Ensembl_mapping_tables)) {
        ## Add mapping to ensembl chromosomes from UCSC tables.
        Ensembl_seqlevel <- .fetch_Ensembl_seqlevels_from_UCSC(
                                           genome,
                                           Ensembl_mapping_tables,
                                           ans$UCSC_seqlevel,
                                           duplicated_contigs,
                                           goldenPath_url=goldenPath_url)
        if (!all(is.na(Ensembl_seqlevel)))
            ans <- cbind(ans, Ensembl_seqlevel=Ensembl_seqlevel,
                              stringsAsFactors=FALSE)
    }
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
        assembly_accession="GCF_000001405.26",
        special_mappings=c(chrM="MT"),
        Ensembl_mapping_tables=c("chromAlias", "ucscToEnsembl"),
        ## The chromInfo table at UCSC contains sequences that belong to
        ## GRCh38.p12 but not to GRCh38. Because we want to map hg38 to
        ## GRCh38 and not to GRCh38.p12, we drop these sequences.
        drop_unmapped=TRUE
    ),

    hg19=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000001405.13",
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
        unmapped_seqs=list(`assembled-molecule`="chrM"),
        Ensembl_mapping_tables=c("chromAlias", "ucscToEnsembl")
    ),

    hg18=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000001405.12",
        special_mappings=c(chr6_cox_hap1="Hs6_111610_36",
                           chr22_h2_hap1="Hs22_111678_36"),
        unmapped_seqs=list(
            `assembled-molecule`="chrM",
            `pseudo-scaffold`=paste0("chr",
                c("5_h2_hap1", "6_qbl_hap2",
                  paste0(c((1:22)[-c(12, 14, 20)], "X"), "_random"))))
        ),

### Chimp
    panTro4=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000001515.5",
        special_mappings=c(chrM="MT"),
        Ensembl_mapping_tables=c("chromAlias", "ucscToEnsembl")
    ),

    panTro3=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCA_000001515.3",
        unmapped_seqs=list(`assembled-molecule`="chrM")
    ),

    panTro2=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000001515.3",
        special_mappings=c(chrM="MT"),
        unmapped_seqs=list(
            `pseudo-scaffold`=c("chr6_hla_hap1", paste0("chr",
                c(1, "2a", "2b", 3:20, 22, "X", "Y"), "_random"),
                "chrUn"))
    ),

### Cow
    bosTau8=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000003055.5",
        special_mappings=c(chrM="MT"),
        Ensembl_mapping_tables="chromAlias"
    ),

    bosTau7=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000003205.5",
        unmapped_seqs=list(`assembled-molecule`="chrM")
    ),

    bosTau6=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000003055.4",
        special_mappings=c(chrM="MT"),
        Ensembl_mapping_tables=c("chromAlias", "ucscToEnsembl")
    ),

### Dog
    canFam3=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000002285.3",
        special_mappings=c(chr1="chr01", chr2="chr02", chr3="chr03",
                           chr4="chr04", chr5="chr05", chr6="chr06",
                           chr7="chr07", chr8="chr08", chr9="chr09",
                           chrM="MT"),
        Ensembl_mapping_tables="chromAlias"
    ),

    canFam2=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM"
    ),

    canFam1=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM"
    ),

### Ferret
    musFur1=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        assembly_accession="GCF_000215625.1",
        Ensembl_mapping_tables="chromAlias"
    ),

### Mouse
    mm10=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000001635.20",
        AssemblyUnits=c("C57BL/6J", "non-nuclear"),
        special_mappings=c(chrM="MT"),
        Ensembl_mapping_tables=c("chromAlias", "ucscToEnsembl")
    ),

    mm9=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000001635.18",
        AssemblyUnits=c("C57BL/6J", "non-nuclear"),
        special_mappings=c(chrM="MT"),
        unmapped_seqs=list(
            `pseudo-scaffold`=paste0("chr",
                c(1, 3:5, 7:9, 13, 16:17, "X", "Y", "Un"), "_random"))
    ),

    mm8=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000001635.15",
        AssemblyUnits="C57BL/6J",
        unmapped_seqs=list(
            `assembled-molecule`="chrM",
            `pseudo-scaffold`=paste0("chr",
                c(1, 5, 7:10, 13, 15, 17, "X", "Y", "Un"), "_random"))
    ),

### Pig
    susScr3=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000003025.5",
        special_mappings=c(chrM="MT"),
        Ensembl_mapping_tables="chromAlias"
    ),

    susScr2=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000003025.3",
        unmapped_seqs=list(`assembled-molecule`="chrM")
    ),

### Rat
    rn6=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000001895.5",
        special_mappings=c(chrM="MT"),
        Ensembl_mapping_tables="chromAlias"
        ## https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/makeDb/doc/rn6.txt
        ## LINES 856:864
        , duplicated_contigs=c("AABR07023006.1", "AABR07022518.1")
    ),

### Rhesus
    rheMac3=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCA_000230795.1",
        unmapped_seqs=list(`assembled-molecule`="chrM")
    ),

    rheMac2=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        assembly_accession="GCF_000002255.3",
        unmapped_seqs=list(`pseudo-scaffold`="chrUr")
    ),

### Chicken
    galGal4=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000002315.3",
        special_mappings=c(chrLGE22C19W28_E50C23="ChrE22C19W28_E50C23",
                           chrLGE64="ChrE64",
                           chrM="MT"),
        Ensembl_mapping_tables=c("chromAlias", "ucscToEnsembl")
    ),

    galGal3=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000002315.2",
        special_mappings=c(chrE22C19W28_E50C23="LGE22C19W28_E50C23",
                           chrE64="LGE64",
                           chrM="MT"),
        unmapped_seqs=list(
            `pseudo-scaffold`=paste0("chr",
                c(1:2, 4:8, 10:13, 16:18, 20, 22, 25, 28, "W", "Z",
                  "E22C19W28_E50C23", "E64", "Un"), "_random"))
    ),

### Stickleback
    gasAcu1=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        #assembly_accession="GCA_000180675.1"
        Ensembl_mapping_tables="chromAlias"
    ),

### Zebra finch
# Check/test this before commenting out. Once everything is ok, update man
# page and commit.
#    taeGut2=list(
#        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
#        circ_seqs="chrM",
#        assembly_accession="GCF_000151805.1",
#        special_mappings=c(chrM="MT")
#    ),
#
#    taeGut1=list(
#        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
#        circ_seqs="chrM"
#    ),

### Zebrafish
    danRer7=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000002035.4",
        special_mappings=c(chrM="MT"),
        Ensembl_mapping_tables="chromAlias"
    ),

    #Too messy!
    #danRer6=list(
    #    FUN="standard_fetch_extended_ChromInfo_from_UCSC",
    #    circ_seqs="chrM",
    #    assembly_accession="GCF_000002035.3",
    #    special_mappings=c(chrM="MT")
    #),

### A. mellifera
    apiMel2=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        assembly_accession="GCF_000002195.1",
        special_mappings=setNames(paste0("LG", 1:16), paste0("Group", 1:16)),
        unmapped_seqs=list(`pseudo-scaffold`="GroupUn")
    ),

### D. melanogaster
    dm6=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000001215.4",
        ## special_mappings=c(chrM="MT"),
        special_mappings=c(chrM="mitochondrion_genome"),
        Ensembl_mapping_tables="chromAlias"
    ),

    dm3=list(
        FUN="standard_fetch_extended_ChromInfo_from_UCSC",
        circ_seqs="chrM",
        assembly_accession="GCF_000001215.2",
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
        assembly_accession="GCF_000002985.1",
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
        assembly_accession="GCF_000146045.2",
        special_mappings=c(chrM="MT"),
        Ensembl_mapping_tables="chromAlias"
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
    .Deprecated("getChromInfoFromUCSC")

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
        assembly_accession=supported_genome$assembly_accession,
        AssemblyUnits=supported_genome$AssemblyUnits,
        special_mappings=supported_genome$special_mappings,
        unmapped_seqs=supported_genome$unmapped_seqs,
        drop_unmapped=supported_genome$drop_unmapped,
        Ensembl_mapping_tables=supported_genome$Ensembl_mapping_tables,
        duplicated_contigs=supported_genome$duplicated_contigs,
        goldenPath_url=goldenPath_url,
        quiet=quiet)
}

