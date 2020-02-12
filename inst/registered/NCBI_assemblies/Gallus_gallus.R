ORGANISM <- "Gallus gallus"

### List of assemblies by date.
### Yep, different genome assemblies can have the same name! (don't ask me why)
### Lookup by genome name will pick-up the first in the list.
ASSEMBLIES <- list(
    list(assembly="Gallus_gallus-2.1",
         date="2006/11/01",
         assembly_accession="GCF_000002315.2",  # galGal3
         circ_seqs="MT"),

    list(assembly="Gallus_gallus-4.0",
         date="2011/11/22",
         assembly_accession="GCA_000002315.2",  # galGal4
         circ_seqs="MT"),

    list(assembly="Gallus_gallus-5.0",
         date="2015/12/16",
         assembly_accession="GCF_000002315.4",  # galGal5
         circ_seqs="MT"),

    list(assembly="GRCg6a",
         date="2018/03/27",
         assembly_accession="GCA_000002315.5",  # galGal6
         circ_seqs="MT"),

    ## Same as above except for the MT sequence!
    list(assembly="GRCg6a",
         date="2018/03/27",
         assembly_accession="GCF_000002315.6",
         circ_seqs="MT")
)

