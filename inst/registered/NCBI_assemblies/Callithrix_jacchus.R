ORGANISM <- "Callithrix jacchus"

### List of assemblies by date.
ASSEMBLIES <- list(
    list(assembly="Callithrix jacchus-3.2",
         date="2010/01/22",
         assembly_accession="GCA_000004665.1",  # calJac3
         circ_seqs=character(0)),

    ## This assembly report has 1992178 contigs! (file to download is 161M)
    list(assembly="ASM83236v1",
         date="2015/01/29",
         assembly_accession="GCA_000832365.1",
         circ_seqs=character(0)),

    ## Note that the SequenceName field is set to 'na' for all the sequences in
    ## the assembly report but getChromInfoFromNCBI("CIEA01") will automatically
    ## set this field to the GenBank accession.
    list(assembly="CIEA01",
         date="2015/08/07",
         assembly_accession="GCA_001269965.1",
         circ_seqs=character(0)),

    list(assembly="ASM275486v1",
         date="2017/11/06",
         assembly_accession="GCA_002754865.1",
         circ_seqs=character(0)),

    list(assembly="Callithrix_jacchus_cj1700_1.0",
         date="2019/11/15",
         assembly_accession="GCA_009663435.1",
         circ_seqs=character(0)),

    list(assembly="Callithrix_jacchus_cj1700_1.1",
         date="2020/05/22",
         assembly_accession="GCF_009663435.1",  # calJac4
         circ_seqs="MT"),

    list(assembly="CJA1912RKC",
         date="2020/06/03",
         assembly_accession="GCA_013373975.1",
         circ_seqs=character(0))
)

