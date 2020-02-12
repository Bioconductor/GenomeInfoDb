GENOME <- "ce10"
ORGANISM <- "Caenorhabditis elegans"
ASSEMBLED_MOLECULES <- paste0("chr", c(as.character(as.roman(1:5)), "X", "M"))
CIRC_SEQS <- "chrM"

### According to UCSC:
###     https://genome.ucsc.edu/cgi-bin/hgGateway?db=ce10
### ce10 "is based on sequence version WS220 deposited into WormBase as of
### October 2010". Problem is that the WS220 assembly doesn't seem to be
### on NCBI.
### OTOH, according to NCBI:
###     https://www.ncbi.nlm.nih.gov/assembly/GCF_000002985.5/
### "ce10" is a synomym for assembly WBcel215 and WBcel215 sequence lengths
### actually match ce10 sequence lengths. So we link ce10 to WBcel215.
NCBI_LINKER <- list(
    assembly_accession="GCF_000002985.5"  # WBcel215
)

