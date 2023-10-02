GENOME <- "loxAfr3"
ORGANISM <- "Loxodonta africana"
ASSEMBLED_MOLECULES <- "chrM"
CIRC_SEQS <- "chrM"

.order_seqlevels <- function(seqlevels)
{
    scaffolds <- paste0("scaffold_", 0:2351)
    ordered_seqlevels <- c(ASSEMBLED_MOLECULES, scaffolds)
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
    assembly_accession="GCF_000001905.1",
    special_mappings=c(chrM="MT"),
    ## For some mysterious reason NCBI and UCSC disagree on the length
    ## of the scaffold_750 sequence. According to the former (Loxafr3.0
    ## assembly) it's 16799 but for the latter (loxAfr3 genome) it's 22731.
    ## This means that the scaffold_750 sequence in one is not the same as
    ## the scaffold_750 sequence in the other so we must not map them.
    ## We also assume that scaffold_750 in loxAfr3 is an "unplaced-scaffold"
    ## (and not an "unlocalized-scaffold" like reported in Loxafr3.0) so it's
    ## no different from all the other scaffolds in the assembly.
    unmapped_seqs=list(`unplaced-scaffold`="scaffold_750")
)

