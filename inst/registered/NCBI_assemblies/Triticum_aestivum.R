ORGANISM <- "Triticum aestivum"

### List of assemblies first by cultivar then by date.
ASSEMBLIES <- list(
    ## ---------- cultivar: Canthatch ---------------

    ## 127921 sequences.
    list(assembly="Can7DL_1.0",
         assembly_level="Contig",
         date="2017/11/16",
         extra_info=c(cultivar="Canthatch K"),
         assembly_accession="GCA_002780475.1",
         circ_seqs=character(0)),

    ## 749802 sequences.
    list(assembly="Canthatch-K",
         assembly_level="Scaffold",
         date="2018/12/03",
         extra_info=c(cultivar="Canthatch"),
         assembly_accession="GCA_900235935.1",
         circ_seqs=character(0)),

    ## 767884 sequences.
    list(assembly="Canthatch-W",
         assembly_level="Scaffold",
         date="2018/12/03",
         extra_info=c(cultivar="Canthatch"),
         assembly_accession="GCA_900235945.1",
         circ_seqs=character(0)),

    ## ---------- cultivar: CDC Landmark ------------

    ## 106140 sequences.
    list(assembly="10wheat_assembly_landmark1",
         assembly_level="Chromosome",
         date="2020/08/22",
         extra_info=c(cultivar="CDC Landmark"),
         assembly_accession="GCA_903995565.1",
         circ_seqs=character(0)),

    ## ---------- cultivar: CDC Stanley -------------

    ## 104381 sequences.
    list(assembly="10wheat_assembly_stanley",
         assembly_level="Chromosome",
         date="2020/08/20",
         extra_info=c(cultivar="CDC Stanley"),
         assembly_accession="GCA_903994155.1",
         circ_seqs=character(0)),

    ## ---------- cultivar: Chinese Spring ----------

    ## 1 sequence.
    list(assembly="ASM21033v1",
         assembly_level="Chromosome",
         date="2010/07/15",
         extra_info=c(cultivar="Chinese Spring"),
         assembly_accession="GCA_000210335.1",
         circ_seqs=character(0)),

    ## 22 sequences.
    list(assembly="iwgsc_refseqv1.0",
         assembly_level="Chromosome",
         date="2018/08/20",
         extra_info=c(cultivar="Chinese Spring (IWGSC RefSeq v1.0)"),
         assembly_accession="GCA_900519105.1",
         circ_seqs=character(0)),

    ## 70316 sequences.
    list(assembly="Triticum_4.0",
         assembly_level="Chromosome",
         date="2020/03/24",
         extra_info=c(cultivar="Chinese Spring"),
         assembly_accession="GCA_002220415.3",
         ## According to the full sequence report the "mitochondrion" sequence
         ## is an unlocalized scaffold. And according to its GenBank page,
         ## the sequence is linear:
         ##   https://www.ncbi.nlm.nih.gov/nuccore/NMPL03070316.1
         circ_seqs="chloroplast"),

    ## 91589 sequences.
    list(assembly="IWGSC CS RefSeq v2.1",
         assembly_level="Chromosome",
         date="2021/05/06",
         extra_info=c(cultivar="Chinese Spring"),
         ## GCA_018294505.1 and GCF_018294505.1 are identical.
         assembly_accession="GCF_018294505.1",
         circ_seqs="MT"),

    ## ---------- cultivar: Fielder -----------------

    ## 3795 sequences.
    list(assembly="wheat_cv_fielder_v1_assembly",
         assembly_level="Chromosome",
         date="2021/07/14",
         extra_info=c(cultivar="Fielder"),
         assembly_accession="GCA_907166925.1",
         circ_seqs=character(0)),

    ## ---------- cultivar: Norin 61 ----------------

    ## 21465 sequences.
    list(assembly="10wheat_assembly_norin61",
         assembly_level="Chromosome",
         date="2020/09/03",
         extra_info=c(cultivar="Norin 61"),
         assembly_accession="GCA_904066035.1",
         circ_seqs=character(0)),

    ## -------------- unknown cultivar --------------

    ## 1 sequence.
    list(assembly="CH Campala Lr22a Pseudomolecule v5",
         assembly_level="Chromosome",
         date="2018/07/17",
         assembly_accession="GCA_900411305.1",
         circ_seqs=character(0)),

    ## 22 sequences.
    list(assembly="Tae_Kariega_v1",
         assembly_level="Chromosome",
         date="2021/12/19",
         assembly_accession="GCA_910594105.1",
         circ_seqs=character(0)),

    ## Something's not working on NCBI's side! (May 9, 2022)
    # > getChromInfoFromNCBI("GCA_920937835.1")
    # Error in function (type, msg, asError = TRUE)  : 
    #   Server denied you to change to the given directory
    list(assembly="Triticum_aestivum_Renan",
         assembly_level="Chromosome",
         date="2022/05/03",
         assembly_accession="GCA_920937835.1",
         circ_seqs=character(0)),

    ## 21 sequences.
    list(assembly="Triticum_aestivum_Renan_v2.1",
         assembly_level="Chromosome",
         date="2022/05/03",
         assembly_accession="GCA_937894285.1",
         circ_seqs=character(0))
)

