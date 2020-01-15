ORGANISM <- "Sus scrofa"

### List of assemblies first by breed then by date.
ASSEMBLIES <- list(
    ## breed: Duroc
    list(genome="Sscrofa9.2",
         assembly_accession="GCF_000003025.3",  # susScr2
         infraspecific_name=c(breed="Duroc"),
         date="2010/02/23",
         circ_seqs=character(0)),

    list(genome="Sscrofa10.2",
         assembly_accession="GCA_000003025.4",  # susScr3
         infraspecific_name=c(breed="Duroc"),
         date="2011/09/07",
         circ_seqs="MT"),

    list(genome="Sscrofa11.1",
         assembly_accession="GCF_000003025.6",  # susScr11
         infraspecific_name=c(breed="Duroc"),
         date="2017/02/07",
         circ_seqs="MT")
)

