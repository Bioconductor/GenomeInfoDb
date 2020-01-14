ORGANISM <- "Canis lupus familiaris"

### List of assemblies first by breed then by date.
### Yep, different genome assemblies can have the same name! (don't ask me why)
### Lookup by genome name will pick-up the first in the list.
ASSEMBLIES <- list(
    ## breed: boxer
    list(genome="CanFam2.0",
         assembly_accession="GCF_000002285.1",
         infraspecific_name=c(breed="boxer"),
         date="2005/07/12",
         circ_seqs=character(0)),

    list(genome="CanFam2.0",
         assembly_accession="GCF_000002285.2",  # canFam2
         infraspecific_name=c(breed="boxer"),
         date="2005/07/12",
         circ_seqs="MT"),

    list(genome="CanFam3.1",
         assembly_accession="GCF_000002285.3",  # canFam3
         infraspecific_name=c(breed="boxer"),
         date="2011/11/02",
         circ_seqs="MT")
)

