GENOME <- "galGal3"
ORGANISM <- "Gallus gallus"
ASSEMBLED_MOLECULES <- paste0("chr",
                              c(1:28, 32, "W", "Z",
                                "E22C19W28_E50C23", "E64",
                                "M"))
CIRC_SEQS <- "chrM"

.random_sequences <- paste0("chr",
                            c(1:2, 4:8, 10:13, 16:18, 20, 22, 25, 28, "W", "Z",
                            "E22C19W28_E50C23", "E64", "Un"), "_random")

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
    assembly_accession="GCF_000002315.2",
    special_mappings=c(chrE22C19W28_E50C23="LGE22C19W28_E50C23",
                       chrE64="LGE64",
                       chrM="MT"),
    unmapped_seqs=list(`pseudo-scaffold`=.random_sequences)
)

