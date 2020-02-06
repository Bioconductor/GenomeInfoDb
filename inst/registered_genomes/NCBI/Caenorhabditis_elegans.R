ORGANISM <- "Caenorhabditis elegans"

### List of assemblies by date.
ASSEMBLIES <- list(
    list(genome="WS144",
         date="2005/06/11",
         assembly_accession="GCA_000002985.1",
         circ_seqs=character(0)),

    list(genome="WS190",
         date="2008/06/27",
         assembly_accession="GCF_000002985.1",  # ce6
         circ_seqs="MT"),

    list(genome="WS195",
         date="2008/11/13",
         assembly_accession="GCF_000002985.4",
         circ_seqs="MT"),

    list(genome="WBcel215",
         date="2012/04/13",
         assembly_accession="GCF_000002985.5",  # ce10
         circ_seqs="MT"),

    list(genome="WBcel235",
         date="2013/02/07",
         extra_info=c(strain="Bristol N2"),
         assembly_accession="GCA_000002985.3",  # ce11
         circ_seqs="MT")
)

