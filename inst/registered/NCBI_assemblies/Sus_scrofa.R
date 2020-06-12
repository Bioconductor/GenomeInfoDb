ORGANISM <- "Sus scrofa"

### List of assemblies first by breed then by date.
ASSEMBLIES <- list(
    ## breed: Duroc
    list(assembly="Sscrofa9.2",
         date="2010/02/23",
         extra_info=c(breed="Duroc"),
         assembly_accession="GCF_000003025.3",  # susScr2
         circ_seqs=character(0)),

    list(assembly="Sscrofa10.2",
         date="2011/09/07",
         extra_info=c(breed="Duroc"),
         assembly_accession="GCA_000003025.4",  # susScr3
         circ_seqs="MT"),

    list(assembly="Sscrofa11.1",
         date="2017/02/07",
         extra_info=c(breed="Duroc"),
         ## This is what Ensembl uses for sus_scrofa in release 99.
         assembly_accession="GCA_000003025.6",  # susScr11
         circ_seqs="MT")
)

