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
GENOME <- "ce6"
ORGANISM <- "Caenorhabditis elegans"
ASSEMBLED_MOLECULES <- paste0("chr", c(as.character(as.roman(1:5)), "X", "M"))

CIRC_SEQS <- "chrM"

