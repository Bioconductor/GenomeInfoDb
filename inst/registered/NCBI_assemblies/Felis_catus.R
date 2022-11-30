ORGANISM <- "Felis catus"

### List of assemblies first by WGS Project, then by date.
ASSEMBLIES <- list(
    ## --- WGS Project: AANG04 ---

    ## 4805 sequences.
    list(assembly="Felis_catus_9.0",
         assembly_level="Chromosome",
         date="2017/11/20",
         extra_info=c(breed="Abyssinian", sex="female"),
         assembly_accession="GCF_000181335.3",  # felCat9
         circ_seqs="MT"),

    ## 4507 sequences.
    list(assembly="felCat9.1_X",
         assembly_level="Chromosome",
         date="2021/10/01",
         extra_info=c(breed="Abyssinian", sex="female"),
         assembly_accession="GCA_000181335.5",
         circ_seqs=character(0)),

    ## --- WGS Project: JAFEKA01 ---

    ## 71 sequences.
    list(assembly="F.catus_Fca126_mat1.0",
         assembly_level="Chromosome",
         date="2021/05/13",
         extra_info=c(sex="female"),
         assembly_accession="GCF_018350175.1",
         circ_seqs="MT")
)

