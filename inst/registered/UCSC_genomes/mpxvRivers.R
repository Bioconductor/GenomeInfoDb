GENOME <- "mpxvRivers"
ORGANISM <- "Monkeypox virus"
ASSEMBLED_MOLECULES <- "NC_063383.1"
CIRC_SEQS <- character(0)

library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

FETCH_ORDERED_CHROM_SIZES <-
    function(goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    chrom_sizes <- GenomeInfoDb:::fetch_chrom_sizes_from_UCSC(GENOME,
                                              from="bigZips",
                                              goldenPath.url=goldenPath.url)
    stopifnot(identical(chrom_sizes[ , "chrom"], ASSEMBLED_MOLECULES))
    chrom_sizes
}

