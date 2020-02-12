GENOME <- "ce6"
ORGANISM <- "Caenorhabditis elegans"
ASSEMBLED_MOLECULES <- paste0("chr", c(as.character(as.roman(1:5)), "X", "M"))
CIRC_SEQS <- "chrM"

NCBI_LINKER <- list(
    assembly_accession="GCF_000002985.1",
    special_mappings=c(chrM="MT")
)

