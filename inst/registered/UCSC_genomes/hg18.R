GENOME <- "hg18"
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
    m33 <- match(m3[ , 3L], paste0("hap", 1:2))
    stopifnot(!anyNA(m33))
    oo3 <- order(m31, m33, m3[ , 2L])
    idx3 <- idx3[oo3]

    idx2 <- which(npart == 2L)
    m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
    m21 <- match(m2[ , 1L], ASSEMBLED_MOLECULES)
    stopifnot(!anyNA(m21))
    stopifnot(all(m2[ , 2L] == "random"))
    oo2 <- order(m21)
    idx2 <- idx2[oo2]

    c(idx1, idx3, idx2)
}

FETCH_ORDERED_CHROM_SIZES <-
    function(goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    chrom_sizes <- GenomeInfoDb:::fetch_chrom_sizes_from_UCSC(GENOME,
                                              goldenPath.url=goldenPath.url)
    oo <- .order_seqlevels(chrom_sizes[ , "chrom"])
    S4Vectors:::extract_data_frame_rows(chrom_sizes, oo)
}

NCBI_LINKER <- list(
    assembly_accession="GCF_000001405.12",
    special_mappings=c(chr6_cox_hap1="Hs6_111610_36",
                       chr22_h2_hap1="Hs22_111678_36"),
    unmapped_seqs=list(
        `assembled-molecule`="chrM",
        `alt-scaffold`=c("chr5_h2_hap1", "chr6_qbl_hap2"),
        `pseudo-scaffold`=
            paste0("chr", c((1:22)[-c(12, 14, 20)], "X"), "_random"))
)

