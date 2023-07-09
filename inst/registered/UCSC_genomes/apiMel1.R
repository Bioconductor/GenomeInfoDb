GENOME <- "apiMel1"
ORGANISM <- "Apis mellifera"
ASSEMBLED_MOLECULES <- character(0)
CIRC_SEQS <- character(0)

library(IRanges)       # for CharacterList()
library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    tmp <- CharacterList(strsplit(seqlevels, ".", fixed=TRUE))
    npart <- lengths(tmp)
    stopifnot(all(npart == 2L))

    m2 <- matrix(unlist(tmp), ncol=2L, byrow=TRUE)
    groups <- paste0("Group", c(1:16, "Un"))
    m21 <- match(m2[, 1L], groups)
    stopifnot(!anyNA(m21))
    m22 <- as.integer(m2[, 2L])
    stopifnot(!anyNA(m22))
    order(m21, m22)
}

FETCH_ORDERED_CHROM_SIZES <-
    function(goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    chrom_sizes <- GenomeInfoDb:::fetch_chrom_sizes_from_UCSC(GENOME,
                                              goldenPath.url=goldenPath.url)
    oo <- .order_seqlevels(chrom_sizes[ , "chrom"])
    S4Vectors:::extract_data_frame_rows(chrom_sizes, oo)
}

