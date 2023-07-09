GENOME <- "bosTau9"
ORGANISM <- "Bos taurus"
ASSEMBLED_MOLECULES <- paste0("chr", c(1:29, "X", "M"))
CIRC_SEQS <- "chrM"

library(IRanges)       # for CharacterList()
library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    tmp <- CharacterList(strsplit(seqlevels, "_"))
    npart <- lengths(tmp)
    stopifnot(all(npart %in% c(1L, 3L)))

    idx1 <- which(npart == 1L)
    stopifnot(length(idx1) == length(ASSEMBLED_MOLECULES))
    oo1 <- match(ASSEMBLED_MOLECULES, seqlevels[idx1])
    stopifnot(!anyNA(oo1))
    idx1 <- idx1[oo1]

    idx3 <- which(npart == 3L)
    m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
    stopifnot(all(m3[ , 1L] == "chrUn"))
    stopifnot(all(m3[ , 2L] == "NW"))
    oo3 <- order(m3[ , 3L])
    idx3 <- idx3[oo3]

    c(idx1, idx3)
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
    assembly_accession="GCA_002263795.2"
)

ENSEMBL_LINKER <- "chromAlias"  # chrM not mapped to Ensembl!

