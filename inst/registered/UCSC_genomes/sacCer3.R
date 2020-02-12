GENOME <- "sacCer3"
ORGANISM <- "Saccharomyces cerevisiae"
ASSEMBLED_MOLECULES <- paste0("chr", c(as.character(as.roman(1:16)), "M"))
CIRC_SEQS <- "chrM"

NCBI_LINKER <- list(
    assembly_accession="GCA_000146045.2",
    special_mappings=c(chrM="MT")
)

ENSEMBL_LINKER <- "chromAlias"

