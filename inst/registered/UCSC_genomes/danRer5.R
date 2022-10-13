GENOME <- "danRer5"
ORGANISM <- "Danio rerio"
ASSEMBLED_MOLECULES <- paste0("chr", c(1:25, "M"))
CIRC_SEQS <- "chrM"

library(IRanges)       # for CharacterList()
library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    tmp <- CharacterList(strsplit(seqlevels, "_"))
    npart <- lengths(tmp)
    stopifnot(all(npart <= 2L))

    idx1 <- which(npart == 1L)
    stopifnot(length(idx1) == length(ASSEMBLED_MOLECULES))
    oo1 <- match(ASSEMBLED_MOLECULES, seqlevels[idx1])
    stopifnot(!anyNA(oo1))
    idx1 <- idx1[oo1]

    idx2 <- which(npart == 2L)
    m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
    stopifnot(all(m2[ , 1L] == "Zv7"))

    idx2_scaffold <- which(substr(m2[ , 2L], 1L, 8L) == "scaffold")
    idx2_NA <- which(substr(m2[ , 2L], 1L, 2L) == "NA")
    stopifnot(length(idx2_scaffold) + length(idx2_NA) == length(idx2))

    m22_scaffold <- m2[idx2_scaffold, 2L]
    suffix <- as.integer(substr(m22_scaffold, 9L, nchar(m22_scaffold)))
    stopifnot(!anyNA(suffix))
    oo2_scaffold <- order(suffix)
    idx2_scaffold <- idx2[idx2_scaffold[oo2_scaffold]]

    m22_NA <- m2[idx2_NA, 2L]
    suffix <- as.integer(substr(m22_NA, 3L, nchar(m22_NA)))
    stopifnot(!anyNA(suffix))
    oo2_NA <- order(suffix)
    idx2_NA <- idx2[idx2_NA[oo2_NA]]

    c(idx1, idx2_scaffold, idx2_NA)
}

FETCH_ORDERED_CHROM_SIZES <-
    function(goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    chrom_sizes <- GenomeInfoDb:::fetch_chrom_sizes_from_UCSC(GENOME,
                                              goldenPath.url=goldenPath.url)
    oo <- .order_seqlevels(chrom_sizes[ , "chrom"])
    S4Vectors:::extract_data_frame_rows(chrom_sizes, oo)
}

