ORGANISM <- "Gallus gallus"

### List of assemblies by date.
### Yep, different genome assemblies can have the same name! (don't ask me why)
### Lookup by genome name will pick-up the first in the list.
ASSEMBLIES <- list(
    list(genome="Gallus_gallus-2.1",
         assembly_accession="GCF_000002315.2",  # galGal3
         date="2006/11/01",
         circ_seqs="MT"),

    list(genome="Gallus_gallus-4.0",
         assembly_accession="GCA_000002315.2",  # galGal4
         date="2011/11/22",
         circ_seqs="MT"),

    list(genome="Gallus_gallus-5.0",
         assembly_accession="GCF_000002315.4",  # galGal5
         date="2015/12/16",
         circ_seqs="MT"),

    list(genome="GRCg6a",
         assembly_accession="GCA_000002315.5",  # galGal6
         date="2018/03/27",
         circ_seqs="MT"),

    ## Same as above except for the MT sequence!
    list(genome="GRCg6a",
         assembly_accession="GCF_000002315.6",
         date="2018/03/27",
         circ_seqs="MT")
)

