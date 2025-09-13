ORGANISM <- "Pan troglodytes"

### List of assemblies by date.
ASSEMBLIES <- list(
    list(assembly="Pan_troglodytes-2.1",
         date="2006/03/16",
         assembly_accession="GCF_000001515.3",  # panTro2
         circ_seqs="MT"),

    ## The sequence names in this one are seriously messed up!
    list(assembly="Pan_troglodytes-2.1.3",
         date="2010/11/15",
         assembly_accession="GCA_000001515.3",  # panTro3
         circ_seqs=character(0)),

    ## 24130 sequences in assembly report but with 2 entries for chromosome 21.
    ## The duplicated entry (RefSeqAccn NC_006488.2) was suddenly added in
    ## Sep 2025 for reasons that are unclear.
    ## Anyways, getChromInfoFromNCBI() will drop this unfortunate addition.
    ## The idea is to keep the original entry because it has length 32,799,110
    ## and this is equal to the length of chr21 in UCSC genome panTro4.
    ## See .get_full_NCBI_chrom_info_from_accession() in
    ## R/getChromInfoFromNCBI.R for the details.
    list(assembly="Pan_troglodytes-2.1.4",
         date="2011/03/25",
         assembly_accession="GCF_000001515.5",  # panTro4
         circ_seqs="MT"),

    list(assembly="Pan_tro 3.0",
         date="2016/05/03",
         assembly_accession="GCF_000001515.7",  # panTro5
         circ_seqs="MT"),

    list(assembly="Clint_PTRv2",
         date="2018/01/19",
         assembly_accession="GCA_002880755.3",  # panTro6
         circ_seqs="MT")
)

