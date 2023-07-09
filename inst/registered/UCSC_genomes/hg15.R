GENOME <- "hg15"
ORGANISM <- "Homo sapiens"
ASSEMBLED_MOLECULES <- paste0("chr", c(1:22, "X", "Y", "M"))
CIRC_SEQS <- "chrM"

.random_sequences <- paste0("chr", c(1:4, 6:13, 15:17, 19, "X", "Y", "Un"),
                            "_random")

library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    ordered_seqlevels <- c(ASSEMBLED_MOLECULES, .random_sequences)
    stopifnot(length(seqlevels) == length(ordered_seqlevels))
    idx <- match(ordered_seqlevels, seqlevels)
    stopifnot(!anyNA(idx))
    idx
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
    assembly_accession="GCF_000001405.8",
    unmapped_seqs=list(
        `assembled-molecule`="chrM",
        `pseudo-scaffold`=.random_sequences)
)

