ORGANISM <- "Macaca fascicularis"

### List of assemblies by date.
ASSEMBLIES <- list(
    list(assembly="Macaca_fascicularis_5.0",
         assembly_level="Chromosome",
         date="2013/06/12",
         extra_info=c(sex="female"),
         assembly_accession="GCA_000364345.1",  # macFas5
         circ_seqs="MT"),

    list(assembly="Macaca_fascicularis_6.0",
         assembly_level="Chromosome",
         date="2020/03/10",
         extra_info=c(sex="male"),
         assembly_accession="GCA_011100615.1",
         circ_seqs=character(0)),

    list(assembly="MFA1912RKS",
         assembly_level="Scaffold",
         date="2020/04/08",
         assembly_accession="GCA_012559485.1",
         circ_seqs=character(0)),

    list(assembly="mfascicularis_v7_1",
         assembly_level="Scaffold",  # seems wrong (contains contigs)
         date="2020/06/05",
         assembly_accession="GCA_903231565.1",
         circ_seqs=character(0)),

    list(assembly="mfascicularis_v7_phased",
         assembly_level="Scaffold",  # seems wrong (contains contigs)
         date="2020/06/05",
         assembly_accession="GCA_903645265.1",
         circ_seqs=character(0)),

    list(assembly="mfascicularis_v7_phased1",
         assembly_level="Scaffold",  # seems wrong (contains contigs)
         date="2020/06/05",
         assembly_accession="GCA_903645275.1",
         circ_seqs=character(0))
)

