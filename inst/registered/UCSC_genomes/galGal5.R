GENOME <- "galGal5"
ORGANISM <- "Gallus gallus"
ASSEMBLED_MOLECULES <- paste0("chr", c(1:28, 30:33, "W", "Z", "LGE64", "M"))
CIRC_SEQS <- "chrM"

library(IRanges)       # for CharacterList()
library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    tmp <- CharacterList(strsplit(seqlevels, "_"))
    npart <- lengths(tmp)
    stopifnot(all(npart %in% c(1L, 3L, 4L)))

    idx1 <- which(npart == 1L)
    stopifnot(length(idx1) == length(ASSEMBLED_MOLECULES))
    oo1 <- match(ASSEMBLED_MOLECULES, seqlevels[idx1])
    stopifnot(!anyNA(oo1))
    idx1 <- idx1[oo1]

    idx4 <- which(npart == 4L)
    m4 <- matrix(unlist(tmp[idx4]), ncol=4L, byrow=TRUE)
    m41 <- match(m4[ , 1L], ASSEMBLED_MOLECULES)
    stopifnot(!anyNA(m41))
    stopifnot(all(m4[ , 2L] == "NT"))
    stopifnot(all(m4[ , 4L] == "random"))
    oo4 <- order(m41, m4[ , 3L])
    idx4 <- idx4[oo4]

    idx3 <- which(npart == 3L)
    m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
    stopifnot(all(m3[ , 1L] == "chrUn"))
    stopifnot(all(m3[ , 2L] == "NT"))
    oo3 <- order(m3[ , 3L])
    idx3 <- idx3[oo3]

    c(idx1, idx4, idx3)
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
    assembly_accession="GCF_000002315.4"
)

ENSEMBL_LINKER <- "chromAlias"

