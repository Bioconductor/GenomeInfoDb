ORGANISM <- "Equus caballus"

ASSEMBLIES <- list(
    ## 9636 sequences.
    list(assembly="EquCab2.0",
         assembly_level="Chromosome",
         date="2007/10/29",
         extra_info=c(breed="thoroughbred"),
         ## GCA_000002305.1 != GCF_000002305.2 (MT only in GCF_000002305.2)
         assembly_accession="GCF_000002305.2",  # equCab2
         circ_seqs="MT"),

    ## 4701 sequences.
    list(assembly="EquCab3.0",
         assembly_level="Chromosome",
         date="2018/01/05",
         extra_info=c(breed="thoroughbred"),
         ## GCA_002863925.1 and GCF_002863925.1 are identical.
         assembly_accession="GCF_002863925.1",  # equCab3
         circ_seqs="MT")
)

