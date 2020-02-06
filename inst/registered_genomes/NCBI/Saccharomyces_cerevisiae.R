ORGANISM <- "Saccharomyces cerevisiae"

### List of assemblies first by strain then by date.
ASSEMBLIES <- list(
    ## strain: BY4742
    list(genome="ASM76643v2",
         date="2014/10/12",
         extra_info=c(strain="BY4742"),
         assembly_accession="GCA_000766435.2",
         circ_seqs=character(0)),

    list(genome="ASM308665v1",
         date="2018/05/03",
         extra_info=c(strain="BY4742"),
         assembly_accession="GCA_003086655.1",
         circ_seqs=character(0)),

    ## strain: S288C
    list(genome="R64",
         date="2014/12/17",
         extra_info=c(strain="S288C"),
         assembly_accession="GCF_000146045.2",  # sacCer3
         circ_seqs="MT"),

    list(genome="ASM205763v1",
         date="2017/03/21",
         extra_info=c(strain="S288C"),
         assembly_accession="GCA_002057635.1",
         circ_seqs=character(0)),

    list(genome="S288c",
         date="2019/07/31",
         extra_info=c(strain="S288C"),
         assembly_accession="GCA_902192305.1",
         circ_seqs=character(0)),

    ## strain: ySR127
    list(genome="ASM105121v1",
         date="2015/07/09",
         extra_info=c(strain="ySR127"),
         assembly_accession="GCA_001051215.1",
         circ_seqs="MT"),

    ## strain: ySR128
    list(genome="ASM432846v1",
         date="2019/03/06",
         extra_info=c(strain="ySR128"),
         assembly_accession="GCA_004328465.1",
         circ_seqs="MT")
)

