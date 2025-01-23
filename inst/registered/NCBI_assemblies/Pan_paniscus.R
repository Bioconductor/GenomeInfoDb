ORGANISM <- "Pan paniscus"

### List of assemblies by date.
ASSEMBLIES <- list(
    list(assembly="panpan1",
         date="2012/04/24",
         extra_info=c(sex="female"),
         assembly_accession="GCA_000258655.1",  # panPan1
         circ_seqs="MT"),

    ## 10275 sequences in assembly report but getChromInfoFromNCBI() will drop
    ## 1 sequence: the duplicated MT sequence (GenBankAccn KT153251.1) that
    ## someone added to this assembly in January 2025. The idea is to keep
    ## the MT sequence with length 16563 because it matches the length of
    ## chrM in UCSC genome panPan2. See .get_NCBI_chrom_info_from_accession()
    ## in R/getChromInfoFromNCBI.R for the details.
    list(assembly="panpan1.1",
         date="2015/08/18",
         extra_info=c(sex="female"),
         assembly_accession="GCF_000258655.2",  # panPan2
         circ_seqs="MT"),

    list(assembly="Mhudiblu_PPA_v0",
         date="2020/05/15",
         extra_info=c(sex="female"),
         assembly_accession="GCF_013052645.1",  # panPan3
         circ_seqs="MT")
)

