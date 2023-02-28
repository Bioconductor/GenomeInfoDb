GENOME <- "hg38"
ORGANISM <- "Homo sapiens"
ASSEMBLED_MOLECULES <- paste0("chr", c(1:22, "X", "Y", "M"))
CIRC_SEQS <- "chrM"

library(IRanges)       # for CharacterList()
library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    tmp <- CharacterList(strsplit(seqlevels, "_"))
    npart <- lengths(tmp)
    stopifnot(all(npart <= 3L))

    idx1 <- which(npart == 1L)
    stopifnot(length(idx1) == length(ASSEMBLED_MOLECULES))
    oo1 <- match(ASSEMBLED_MOLECULES, seqlevels[idx1])
    stopifnot(!anyNA(oo1))
    idx1 <- idx1[oo1]

    idx3 <- which(npart == 3L)
    m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
    m31 <- match(m3[ , 1L], ASSEMBLED_MOLECULES)
    stopifnot(!anyNA(m31))
    m33 <- match(m3[ , 3L], c("alt", "random", "fix"))
    stopifnot(!anyNA(m33))
    GenBankAccn_prefixes <- c("GL",  # alt-scaffold or unlocalized-scaffold
                              "JH",  # alt-scaffold
                              "KB",  # alt-scaffold
                              "KI",  # alt-scaffold or unlocalized-scaffold
                              "KN",  # fix-patch or novel-patch
                              "KQ",  # fix-patch or novel-patch
                              "KV",  # fix-patch or novel-patch
                              "KZ",  # fix-patch or novel-patch
                              "ML",  # fix-patch or novel-patch
                              "MU")  # fix-patch or novel-patch
    m32 <- match(substr(m3[ , 2L], 1L, 2L), GenBankAccn_prefixes)
    stopifnot(!anyNA(m32))
    is_alt_scaffold <- m33 == 1L & m32 <= 4L
    is_unlocalized_scaffold <- m33 == 2L & m32 <= 4L
    is_fix_patch <- m33 == 3L & m32 >= 5L
    is_novel_patch <- m33 == 1L & m32 >= 5L
    sequence_role <- integer(length(idx3))
    sequence_role[is_alt_scaffold] <- 1L
    sequence_role[is_unlocalized_scaffold] <- 2L
    sequence_role[is_fix_patch] <- 3L
    sequence_role[is_novel_patch] <- 4L
    stopifnot(all(sequence_role != 0L))
    oo3 <- order(sequence_role, m31, m3[ , 2L])
    idx3 <- idx3[oo3]
    npatch <- sum(is_fix_patch) + sum(is_novel_patch)
    idx3_patch <- tail(idx3, n=npatch)
    if (npatch != 0L)
        idx3 <- head(idx3, n=-npatch)

    idx2 <- which(npart == 2L)
    m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
    stopifnot(all(m2[ , 1L] == "chrUn"))
    oo2 <- order(m2[ , 2L])
    idx2 <- idx2[oo2]

    c(idx1, idx3, idx2, idx3_patch)
}

FETCH_ORDERED_CHROM_SIZES <-
    function(goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    chrom_sizes <- GenomeInfoDb:::fetch_chrom_sizes_from_UCSC(GENOME,
                                              goldenPath.url=goldenPath.url)
    oo <- .order_seqlevels(chrom_sizes[ , "chrom"])
    S4Vectors:::extract_data_frame_rows(chrom_sizes, oo)
}

#NCBI_LINKER <- list(
#    assembly_accession="GCF_000001405.26",  # GRCh38
#    ## The chromInfo table at UCSC contains sequences that belong to
#    ## GRCh38.p12 but not to GRCh38. Because we want to map hg38 to
#    ## GRCh38 and not to GRCh38.p12, we drop these sequences.
#    drop_unmapped=TRUE
#)

NCBI_LINKER <- list(
    assembly_accession="GCF_000001405.40",  # GRCh38.p14
    ## Note that FIX PATCH sequences HG107_PATCH and HG1311_PATCH
    ## in GRCh38.p13 have been replaced with FIX PATCH sequences
    ## HG107_HG2565_PATCH and HG1311_HG2539_PATCH in GRCh38.p14!
    ## However, when UCSC sneakily modified their hg38 genome towards the
    ## end of January 2023 to base it on GRCh38.p14 instead of GRCh38.p13,
    ## they did the following:
    ## - They kept the old HG107_PATCH and HG1311_PATCH sequences, even
    ##   though these sequences do not belong to GRCh38.p14. These are
    ##   named chr11_KQ759759v1_fix and chr22_KQ759762v1_fix in hg38.
    ## - They also included their replacements, HG107_HG2565_PATCH and
    ##   HG1311_HG2539_PATCH. These are named chr11_KQ759759v2_fix and
    ##   chr22_KQ759762v2_fix in hg38.
    ## This means that hg38 has 2 extra sequences w.r.t. GRCh38.p14 (711
    ## sequences in hg38 vs 709 sequences in GRCh38.p14).
    unmapped_seqs=list(`fix-patch`=c("chr11_KQ759759v1_fix",
                                     "chr22_KQ759762v1_fix"))
)

### Sequences not in the original GRCh38 are not mapped to Ensembl!
ENSEMBL_LINKER <- "ucscToEnsembl"

