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
GENOME <- "dm3"
ORGANISM <- "Drosophila melanogaster"
ASSEMBLED_MOLECULES <- paste0("chr", c("2L", "2R", "3L", "3R", "4", "X", "M"))

CIRC_SEQS <- "chrM"

library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    het <- c("chr2LHet", "chr2RHet", "chr3LHet", "chr3RHet",
             "chrXHet", "chrYHet")
    ordered_seqlevels <- c(ASSEMBLED_MOLECULES, het, "chrU", "chrUextra")
    stopifnot(length(seqlevels) == length(ordered_seqlevels))
    idx <- match(ordered_seqlevels, seqlevels)
    stopifnot(!anyNA(idx))
    idx
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
    assembly_accession="GCF_000001215.2",
    special_mappings=c(chrM="MT", chrU="Un"),
    unmapped_seqs=list(`pseudo-scaffold`="chrUextra")
)

