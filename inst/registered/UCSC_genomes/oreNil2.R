GENOME <- "oreNil2"
ORGANISM <- "Oreochromis niloticus"
ASSEMBLED_MOLECULES <- paste0("chr",
    c(paste0("LG", c(1:7, "8-24", 9:15, "16-21", 17:20, 22:23)), "M"))
CIRC_SEQS <- "chrM"

.order_seqlevels <- function(seqlevels)
{
    idx1 <- which(seqlevels %in% ASSEMBLED_MOLECULES)
    stopifnot(length(idx1) == length(ASSEMBLED_MOLECULES))
    oo1 <- match(ASSEMBLED_MOLECULES, seqlevels[idx1])
    stopifnot(!anyNA(oo1))
    idx1 <- idx1[oo1]

    idx2 <- which(substr(seqlevels, 1, 8) == "AERX0107")
    idx3 <- which(substr(seqlevels, 1, 4) == "GL83")
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
    assembly_accession="GCA_000188235.2",
    special_mappings=c(chrM="MT")
)

ENSEMBL_LINKER <- "chromAlias"

