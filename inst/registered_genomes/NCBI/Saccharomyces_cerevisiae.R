ORGANISM <- "Saccharomyces cerevisiae"

### List of assemblies first by strain then by date.
ASSEMBLIES <- list(
    ## strain: BY4742
    list(genome="ASM76643v2",
         assembly_accession="GCA_000766435.2",
         infraspecific_name=c(strain="BY4742"),
         date="2014/10/12",
         circ_seqs=character(0)),

    list(genome="ASM308665v1",
         assembly_accession="GCA_003086655.1",
         infraspecific_name=c(strain="BY4742"),
         date="2018/05/03",
         circ_seqs=character(0)),

    ## strain: S288C
    list(genome="R64",
         assembly_accession="GCF_000146045.2",  # sacCer3
         infraspecific_name=c(strain="S288C"),
         date="2014/12/17",
         circ_seqs="MT"),

    list(genome="ASM205763v1",
         assembly_accession="GCA_002057635.1",
         infraspecific_name=c(strain="S288C"),
         date="2017/03/21",
         circ_seqs=character(0)),

    list(genome="S288c",
         assembly_accession="GCA_902192305.1",
         infraspecific_name=c(strain="S288C"),
         date="2019/07/31",
         circ_seqs=character(0)),

    ## strain: ySR127
    list(genome="ASM105121v1",
         assembly_accession="GCA_001051215.1",
         infraspecific_name=c(strain="ySR127"),
         date="2015/07/09",
         circ_seqs="MT"),

    ## strain: ySR128
    list(genome="ASM432846v1",
         assembly_accession="GCA_004328465.1",
         infraspecific_name=c(strain="ySR128"),
         date="2019/03/06",
         circ_seqs="MT")
)

