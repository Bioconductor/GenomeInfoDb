GENOME <- "galGal2"
ORGANISM <- "Gallus gallus"
ASSEMBLED_MOLECULES <- paste0("chr",
                              c(1:24, 26:28, 32, "W", "Z",
                                "E26C13", "E50C23", "E22C19W28", "E64",
                                "M"))
CIRC_SEQS <- "chrM"

library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    random <- paste0("chr", c(1:8, 10:11, 13, 16, 24, 27:28, 32, "W", "Z"),
                     "_random")
    ordered_seqlevels <- c(ASSEMBLED_MOLECULES, random, "chrUn")
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

