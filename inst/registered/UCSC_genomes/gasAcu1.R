GENOME <- "gasAcu1"
ORGANISM <- "Gasterosteus aculeatus"
ASSEMBLED_MOLECULES <- paste0("chr", c(as.character(as.roman(1:21)), "M"))
CIRC_SEQS <- "chrM"

library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    ordered_seqlevels <- c(ASSEMBLED_MOLECULES, "chrUn")
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

### UCSC claims that gasAcu1 is based on the GCA_000180675.1 assembly:
###     https://genome.ucsc.edu/cgi-bin/hgGateway?db=gasAcu1
### but none of the sequences in gasAcu1 actually corresponds to a
### sequence in GCA_000180675.1 (which is made of contigs only).
#NCBI_LINKER <- list(
#    assembly_accession="GCA_000180675.1",
#    unmapped_seqs=list(`assembled-molecule`="chrM",
#                       `pseudo-scaffold`="chrUn")
#)

ENSEMBL_LINKER <- "chromAlias"  # chrM not mapped to Ensembl!

