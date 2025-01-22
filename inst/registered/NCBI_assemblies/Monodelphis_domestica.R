ORGANISM <- "Monodelphis domestica"

### List of assemblies by date.
ASSEMBLIES <- list(
    ## 5138 sequences.
    ## Note that this assembly used to have only 5016 sequences but, on
    ## 2022/06/27, someone decided to add 4022 sequences to it, hence
    ## increasing the number of sequences to 9038. The added sequences are
    ## all the sequences after the MT line in the assembly report at
    ## https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/295/GCF_000002295.2_MonDom5/GCF_000002295.2_MonDom5_assembly_report.txt
    ## Yeah they did this. Instead of submitting a new version of the assembly
    ## like everybody else! Also who would have known that NCBI allows
    ## submitters to silently modify an assembly that they submitted years
    ## before without any mention on the assembly landing page, without
    ## assigning it a new RefSeq assembly accession, and without even bothering
    ## to bump the assembly date!
    ## Wait what? Now on 2025/01/09 they did it again! This time they removed
    ## 3900 sequences from the assembly, reducing the number of sequences to
    ## 5138. And again, this modified assembly is still associated with GenBank
    ## and RefSeq accessions GCA_000002295.1 and GCF_000002295.2, respectively.
    ## And of course, nobody seemed to think that updating the date of the
    ## would be a good thing. Date of assembly is still 2007-01-25 on assembly
    ## report. Great job "The Genome Sequencing Platform, The Genome Assembly
    ## Team"! Yes, that's the name of the submitter of this assembly. That's
    ## how whoever submitted this assembly identifies themselves... SMH
    list(assembly="MonDom5",
         date="2007/01/25",
         assembly_accession="GCF_000002295.2",  # monDom5
         circ_seqs="MT")
)

