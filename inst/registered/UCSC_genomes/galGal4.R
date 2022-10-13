GENOME <- "galGal4"
ORGANISM <- "Gallus gallus"
ASSEMBLED_MOLECULES <- paste0("chr",
                              c(1:28, 32, "W", "Z",
                                "LGE22C19W28_E50C23", "LGE64",
                                "M"))
CIRC_SEQS <- "chrM"

library(IRanges)       # for CharacterList()
library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    idx1 <- match(ASSEMBLED_MOLECULES, seqlevels)
    stopifnot(!anyNA(idx1))

    tmp <- CharacterList(strsplit(seqlevels, "_"))
    npart <- lengths(tmp)
    stopifnot(all(npart <= 3L))

    idx3 <- which(npart == 3L)
    m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
    assembled <- ASSEMBLED_MOLECULES
    assembled[assembled == "chrLGE22C19W28_E50C23"] <- "chrLGE22C19W28"
    m31 <- match(m3[ , 1L], assembled)
    stopifnot(!anyNA(m31))
    stopifnot(all(m3[ , 3L] == "random"))
    oo3 <- order(m31, m3[ , 2L])
    idx3 <- idx3[oo3]

    idx2 <- which(npart == 2L & seqlevels != "chrLGE22C19W28_E50C23")
    m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
    stopifnot(all(m2[ , 1L] == "chrUn"))
    oo2 <- order(m2[ , 2L])
    idx2 <- idx2[oo2]

    c(idx1, idx3, idx2)
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
    assembly_accession="GCA_000002315.2",
    special_mappings=c(chrLGE22C19W28_E50C23="ChrE22C19W28_E50C23",
                       chrLGE64="ChrE64",
                       chrM="MT")
)

ENSEMBL_LINKER <- "ucscToEnsembl"

