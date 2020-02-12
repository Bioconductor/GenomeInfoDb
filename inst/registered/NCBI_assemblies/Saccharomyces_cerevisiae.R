ORGANISM <- "Saccharomyces cerevisiae"

### List of assemblies first by strain then by date.
ASSEMBLIES <- list(
    ## strain: BY4742
    list(assembly="ASM76643v2",
         date="2014/10/12",
         extra_info=c(strain="BY4742"),
         assembly_accession="GCA_000766435.2",
         circ_seqs=character(0)),

    list(assembly="ASM308665v1",
         date="2018/05/03",
         extra_info=c(strain="BY4742"),
         assembly_accession="GCA_003086655.1",
         circ_seqs=character(0)),

    ## strain: S288C
    list(assembly="R64",      # called R64-1-1 in Ensembl
         date="2014/12/17",
         extra_info=c(strain="S288C"),
         ## GCA_000146045.2 and GCF_000146045.2 are equivalent but
         ## we use the former because that's what Ensembl uses for
         ## saccharomyces_cerevisiae in release 99. Strangely they
         ## call this assembly R64-1-1 (like UCSC).
         assembly_accession="GCA_000146045.2",  # sacCer3
         circ_seqs="MT",
         ## MT is called Mito at Ensembl. MT (GenBank=KP263414.1,
         ## RefSeq=NC_001224.1) and Mito (GenBank=AJ011856.1) are
         ## exactly the same DNA sequence (I checked that). Don't
         ## ask me why we need 2 GenBank entries for the same sequence
         ## or why only the latter is associated with a RefSeq accession!
         NCBI2Ensembl_special_mappings=c(MT="Mito")),

    list(assembly="ASM205763v1",
         date="2017/03/21",
         extra_info=c(strain="S288C"),
         assembly_accession="GCA_002057635.1",
         circ_seqs=character(0)),

    list(assembly="S288c",
         date="2019/07/31",
         extra_info=c(strain="S288C"),
         assembly_accession="GCA_902192305.1",
         circ_seqs=character(0)),

    ## strain: ySR127
    list(assembly="ASM105121v1",
         date="2015/07/09",
         extra_info=c(strain="ySR127"),
         assembly_accession="GCA_001051215.1",
         circ_seqs="MT"),

    ## strain: ySR128
    list(assembly="ASM432846v1",
         date="2019/03/06",
         extra_info=c(strain="ySR128"),
         assembly_accession="GCA_004328465.1",
         circ_seqs="MT")
)

