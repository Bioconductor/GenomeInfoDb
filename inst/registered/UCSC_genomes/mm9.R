GENOME <- "mm9"
ORGANISM <- "Mus musculus"
ASSEMBLED_MOLECULES <- paste0("chr", c(1:19, "X", "Y", "M"))
CIRC_SEQS <- "chrM"

library(IRanges)       # for CharacterList()
library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    idx_chrUn <- match("chrUn_random", seqlevels)
    stopifnot(!anyNA(idx_chrUn))

    tmp <- CharacterList(strsplit(seqlevels, "_"))
    npart <- lengths(tmp)
    stopifnot(all(npart <= 2L))

    idx1 <- which(npart == 1L)
    stopifnot(length(idx1) == length(ASSEMBLED_MOLECULES))
    oo1 <- match(ASSEMBLED_MOLECULES, seqlevels[idx1])
    stopifnot(!anyNA(oo1))
    idx1 <- idx1[oo1]

    idx2 <- which(npart == 2L & seqlevels != "chrUn_random")
    m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
    m21 <- match(m2[ , 1L], ASSEMBLED_MOLECULES)
    stopifnot(!anyNA(m21))
    stopifnot(all(m2[ , 2L] == "random"))
    oo2 <- order(m21)
    idx2 <- idx2[oo2]

    c(idx1, idx2, idx_chrUn)
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
    assembly_accession="GCF_000001635.18",
    AssemblyUnits=c("C57BL/6J", "non-nuclear"),
    unmapped_seqs=list(
        `pseudo-scaffold`=
            paste0("chr", c(1, 3:5, 7:9, 13, 16:17, "X", "Y", "Un"), "_random")
    )
)

