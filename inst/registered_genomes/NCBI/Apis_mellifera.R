ORGANISM <- "Apis mellifera"

### List of assemblies first by submitter then by date.
ASSEMBLIES <- list(
    ## submitter: Human Genome Sequencing Center
    list(genome="Amel_2.0",
         date="2005/01/25",
         extra_info=c(strain="DH4"),
         assembly_accession="GCF_000002195.1",  # apiMel2
         circ_seqs=character(0)),

    list(genome="Amel_4.0",
         date="2005/05/05",
         extra_info=c(strain="DH4"),
         assembly_accession="GCF_000002195.2",
         circ_seqs=character(0)),

    list(genome="Amel_4.0",
         date="2005/05/05",
         extra_info=c(strain="DH4"),
         assembly_accession="GCF_000002195.3",
         circ_seqs="MT"),

    list(genome="Amel_4.5",
         date="2011/01/14",
         extra_info=c(strain="DH4"),
         assembly_accession="GCF_000002195.4",
         circ_seqs="MT"),

    ## submitter: INRA
    list(genome="INRA_AMelMel_1.0",
         date="2018/07/11",
         extra_info=c(sex="male", submitter="INRA"),
         assembly_accession="GCA_003314205.1",
         circ_seqs=character(0)),

    ## submitter: Uppsala University
    list(genome="Amel_HAv3",
         date="2018/06/19",
         extra_info=c(strain="DH4", sex="male", submitter="Uppsala University"),
         assembly_accession="GCA_003254395.1",
         circ_seqs="MT"),

    list(genome="Amel_HAv3.1",
         date="2018/09/10",
         extra_info=c(strain="DH4", sex="male", submitter="Uppsala University"),
         assembly_accession="GCF_003254395.2",
         circ_seqs="MT")  # MT length set to NA in assembly report!
)

