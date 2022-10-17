GENOME <- "mpxvRivers"
ORGANISM <- "Monkeypox virus"
ASSEMBLED_MOLECULES <- "NC_063383.1"
CIRC_SEQS <- character(0)

library(GenomeInfoDb)

FETCH_ORDERED_CHROM_SIZES <-
    function(goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    filename <- paste0(GENOME, ".chrom.sizes.txt")
    url <- paste(goldenPath.url, GENOME, "bigZips", filename, sep="/")
    col2class <- c(chrom="character", size="integer")
    ans <- GenomeInfoDb:::fetch_table_from_url(url,
                                      colnames=names(col2class),
                                      col2class=col2class)
    stopifnot(identical(ans[ , "chrom"], ASSEMBLED_MOLECULES))
    ans
}

