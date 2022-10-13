GENOME <- "panTro2"
ORGANISM <- "Pan troglodytes"
ASSEMBLED_MOLECULES <- paste0("chr", c(1, "2a", "2b", 3:22, "X", "Y", "M"))
CIRC_SEQS <- "chrM"

library(IRanges)       # for CharacterList()
library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    idx_chrUn <- match("chrUn", seqlevels)
    stopifnot(!anyNA(idx_chrUn))

    tmp <- CharacterList(strsplit(seqlevels, "_"))
    npart <- lengths(tmp)
    stopifnot(all(npart <= 3L))

    idx1 <- which(npart == 1L & seqlevels != "chrUn")
    stopifnot(length(idx1) == length(ASSEMBLED_MOLECULES))
    oo1 <- match(ASSEMBLED_MOLECULES, seqlevels[idx1])
    stopifnot(!anyNA(oo1))
    idx1 <- idx1[oo1]

    idx3 <- which(npart == 3L)
    stopifnot(length(idx3) == 1L, seqlevels[idx3] == "chr6_hla_hap1")

    idx2 <- which(npart == 2L)
    m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
    m21 <- match(m2[ , 1L], ASSEMBLED_MOLECULES)
    stopifnot(!anyNA(m21))
    stopifnot(all(m2[ , 2L] == "random"))
    oo2 <- order(m21)
    idx2 <- idx2[oo2]

    c(idx1, idx3, idx2, idx_chrUn)
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
    assembly_accession="GCF_000001515.3",
    special_mappings=c(chrM="MT"),
    unmapped_seqs=list(
        `pseudo-scaffold`=
            c("chr6_hla_hap1",
              paste0("chr", c(1, "2a", "2b", 3:20, 22, "X", "Y"), "_random"),
              "chrUn")
    )
)

