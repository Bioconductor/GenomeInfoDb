GENOME <- "taeGut1"
ORGANISM <- "Taeniopygia guttata"
ASSEMBLED_MOLECULES <- paste0("chr",
                              c(1, "1A", "1B", 2:4, "4A", 5:28,
                                "Z", "LGE22", "LG2", "LG5",
                                "M"))
CIRC_SEQS <- "chrM"

library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    random <- paste0(head(ASSEMBLED_MOLECULES, n=33L), "_random")
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

ENSEMBL_LINKER <- "chromAlias"

