GENOME <- "dm3"
ORGANISM <- "Drosophila melanogaster"
ASSEMBLED_MOLECULES <- paste0("chr", c("2L", "2R", "3L", "3R", "4", "X", "M"))
CIRC_SEQS <- "chrM"

library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    het <- c("chr2LHet", "chr2RHet", "chr3LHet", "chr3RHet",
             "chrXHet", "chrYHet")
    ordered_seqlevels <- c(ASSEMBLED_MOLECULES, het, "chrU", "chrUextra")
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
    assembly_accession="GCF_000001215.2",
    special_mappings=c(chrM="MT", chrU="Un"),
    unmapped_seqs=list(`pseudo-scaffold`="chrUextra")
)

