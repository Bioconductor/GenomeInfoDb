### Should work as a standlone, self-contained script.
### Must define at least:
###   o GENOME:              Single non-empty string.
###   o ORGANISM:            Single non-empty string.
###   o ASSEMBLED_MOLECULES: Character vector with no NAs, no empty strings,
###                          and no duplicates.
### Can also define:
###   o CIRC_SEQS:           Character vector (subset of ASSEMBLED_MOLECULES).
###   o GET_CHROM_SIZES:     Function with 1 argument. Must return a 2-column
###                          data.frame with columns "chrom" and "size".
GENOME <- "hg18"
ORGANISM <- "Homo sapiens"
ASSEMBLED_MOLECULES <- paste0("chr", c(1:22, "X", "Y", "M"))

CIRC_SEQS <- "chrM"

library(IRanges)       # for CharacterList()
library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    tmp <- CharacterList(strsplit(seqlevels, "_"))
    npart <- lengths(tmp)
    stopifnot(all(npart <= 3L))

    idx1 <- which(npart == 1L)
    stopifnot(length(idx1) == length(ASSEMBLED_MOLECULES))
    oo1 <- match(ASSEMBLED_MOLECULES, seqlevels[idx1])
    stopifnot(!anyNA(oo1))
    idx1 <- idx1[oo1]

    idx3 <- which(npart == 3L)
    m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
    m31 <- match(m3[ , 1L], ASSEMBLED_MOLECULES)
    stopifnot(!anyNA(m31))
    m33 <- match(m3[ , 3L], paste0("hap", 1:2))
    stopifnot(!anyNA(m33))
    oo3 <- order(m31, m33, m3[ , 2L])
    idx3 <- idx3[oo3]

    idx2 <- which(npart == 2L)
    m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
    m21 <- match(m2[ , 1L], ASSEMBLED_MOLECULES)
    stopifnot(!anyNA(m21))
    stopifnot(all(m2[ , 2L] == "random"))
    oo2 <- order(m21)
    idx2 <- idx2[oo2]

    c(idx1, idx3, idx2)
}

GET_CHROM_SIZES <- function(goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    chrom_sizes <- GenomeInfoDb:::fetch_chrom_sizes_from_UCSC(GENOME,
                                              goldenPath.url=goldenPath.url)
    oo <- .order_seqlevels(chrom_sizes[ , "chrom"])
    S4Vectors:::extract_data_frame_rows(chrom_sizes, oo)
}

