ORGANISM <- "Monodelphis domestica"

### List of assemblies by date.
ASSEMBLIES <- list(
    ## 9038 sequences. Note that this assembly used to have only 5016
    ## sequences but, on 2022/06/27, someone added 4022 sequences to it.
    ## The new sequences are all the sequences after the MT line in
    ## GCF_000002295.2_MonDom5_assembly_report.txt. Yeah they did this,
    ## instead of submitting a new version of the assembly like everybody
    ## else! Also who would have known that NCBI allows submitters to
    ## silently modify an assembly that they submitted years before without
    ## any mention on the assembly landing page, without assigning it a new
    ## RefSeq assembly accession, and without even bumping the assembly date!
    list(assembly="MonDom5",
         date="2007/01/25",
         assembly_accession="GCF_000002295.2",  # monDom5
         circ_seqs="MT")
)

