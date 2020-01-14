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
###   o NCBI_LINKER:         Named list.
GENOME <- "susScr3"
ORGANISM <- "Sus scrofa"
ASSEMBLED_MOLECULES <- paste0("chr", c(1:18, "X", "Y", "M"))

CIRC_SEQS <- "chrM"

library(IRanges)       # for CharacterList()
library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    tmp <- CharacterList(strsplit(seqlevels, "-"))
    npart <- lengths(tmp)
    stopifnot(all(npart <= 2L))

    idx1 <- which(npart == 1L)
    stopifnot(length(idx1) == length(ASSEMBLED_MOLECULES))
    oo1 <- match(ASSEMBLED_MOLECULES, seqlevels[idx1])
    stopifnot(!anyNA(oo1))
    idx1 <- idx1[oo1]

    idx2 <- which(npart == 2L)
    oo2 <- order(seqlevels[idx2])
    idx2 <- idx2[oo2]

    c(idx1, idx2)
}

GET_CHROM_SIZES <- function(goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    chrom_sizes <- GenomeInfoDb:::fetch_chrom_sizes_from_UCSC(GENOME,
                                              goldenPath.url=goldenPath.url)
    oo <- .order_seqlevels(chrom_sizes[ , "chrom"])
    S4Vectors:::extract_data_frame_rows(chrom_sizes, oo)
}

### Valid NCBI_LINKER components:
### - assembly_accession: single non-empty string.
### - AssemblyUnits: character vector.
### - special_mappings: named character vector.
### - unmapped_seqs: named list of character vectors.
### - drop_unmapped: TRUE or FALSE.
NCBI_LINKER <- list(
    assembly_accession="GCF_000003025.5",
    special_mappings=c(chrM="MT")
)

