GENOME <- "mm10"
ORGANISM <- "Mus musculus"
ASSEMBLED_MOLECULES <- paste0("chr", c(1:19, "X", "Y", "M"))
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

    ## fix patches
    fix_patch_idx <- grep("_fix$", seqlevels)
    m3 <- matrix(unlist(tmp[fix_patch_idx]), ncol=3L, byrow=TRUE)
    m31 <- match(m3[ , 1L], ASSEMBLED_MOLECULES)
    stopifnot(!anyNA(m31))
    oo3 <- order(m31, m3[ , 2L])
    fix_patch_idx <- fix_patch_idx[oo3]

    ## novel patches
    novel_patches <- c("chr1_KK082441_alt",
                       "chr11_KZ289073_alt",
                       "chr11_KZ289074_alt",
                       "chr11_KZ289075_alt",
                       "chr11_KZ289077_alt",
                       "chr11_KZ289078_alt",
                       "chr11_KZ289079_alt",
                       "chr11_KZ289080_alt",
                       "chr11_KZ289081_alt")
    novel_patch_idx <- match(novel_patches, seqlevels)
    stopifnot(!anyNA(novel_patch_idx))

    ## alt scaffolds and unlocalized scaffolds
    idx3 <- setdiff(which(npart == 3L), c(fix_patch_idx, novel_patch_idx))
    m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
    m31 <- match(m3[ , 1L], c(ASSEMBLED_MOLECULES, "chrna"))
    stopifnot(!anyNA(m31))
    m33 <- match(m3[ , 3L], c("alt", "random"))
    stopifnot(!anyNA(m33))
    oo3 <- order(m33, m31, m3[ , 2L])
    idx3 <- idx3[oo3]

    ## unplaced scaffolds
    idx2 <- which(npart == 2L)
    m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
    stopifnot(all(m2[ , 1L] == "chrUn"))
    oo2 <- order(m2[ , 2L])
    idx2 <- idx2[oo2]

    c(idx1, idx3, idx2, fix_patch_idx, novel_patch_idx)
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
    assembly_accession="GCF_000001635.26"
)

ENSEMBL_LINKER <- "ucscToEnsembl"

