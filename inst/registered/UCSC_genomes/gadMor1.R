GENOME <- "gadMor1"
ORGANISM <- "Gadus morhua"
ASSEMBLED_MOLECULES <- CIRC_SEQS <- "chrM"

.order_seqlevels <- function(seqlevels)
{
    idx1 <- match(ASSEMBLED_MOLECULES, seqlevels)
    stopifnot(!anyNA(idx1))

    idx2 <- which(substr(seqlevels, 1, 6) == "CAEA01")
    idx3 <- which(substr(seqlevels, 1, 3) == "HE5")
    stopifnot(length(idx1) + length(idx2) + length(idx3) == length(seqlevels))

    oo2 <- order(seqlevels[idx2])
    idx2 <- idx2[oo2]

    oo3 <- order(seqlevels[idx3])
    idx3 <- idx3[oo3]

    c(idx1, idx2, idx3)
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
    assembly_accession="GCA_000231765.1",
    unmapped_seqs=list(`assembled-molecule`="chrM")
)

ENSEMBL_LINKER <- "ucscToEnsembl"

