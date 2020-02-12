ORGANISM <- "Canis lupus familiaris"

### List of assemblies first by breed then by date.
### Yep, different genome assemblies can have the same name! (don't ask me why)
### Lookup by genome name will pick-up the first in the list.
ASSEMBLIES <- list(
    ## breed: boxer
    list(assembly="CanFam2.0",
         date="2005/07/12",
         extra_info=c(breed="boxer"),
         assembly_accession="GCF_000002285.1",
         circ_seqs=character(0)),

    list(assembly="CanFam2.0",
         date="2005/07/12",
         extra_info=c(breed="boxer"),
         assembly_accession="GCF_000002285.2",  # canFam2
         circ_seqs="MT"),

    list(assembly="CanFam3.1",
         date="2011/11/02",
         extra_info=c(breed="boxer"),
         assembly_accession="GCF_000002285.3",  # canFam3
         circ_seqs="MT")
)

