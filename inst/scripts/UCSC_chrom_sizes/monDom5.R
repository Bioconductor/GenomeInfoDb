### Should work as a standlone, self-contained script.
### Must define at least:
###   o ASSEMBLY:            Single non-empty string.
###   o ASSEMBLED_MOLECULES: Character vector with no NAs, no empty strings,
###                          and no duplicates.
### Can also define:
###   o GET_CHROM_SIZES:     Function with 1 argument. Must return a 2-column
###                          data.frame with columns "chrom" and "size".
ASSEMBLY <- "monDom5"
ASSEMBLED_MOLECULES <- paste0("chr", c(1:8, "X", "M", "Un"))

