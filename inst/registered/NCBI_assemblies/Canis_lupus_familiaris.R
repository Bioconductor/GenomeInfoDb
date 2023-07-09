ORGANISM <- "Canis lupus familiaris"

### List of assemblies first by breed then by date.
### Yep, different genome assemblies can have the same name! (don't ask me why)
### Lookup by genome name will pick-up the first in the list.
ASSEMBLIES <- list(
    ## breed: boxer
    list(assembly="CanFam2.0",
         date="2005/07/12",
         extra_info=c(breed="boxer"),
         assembly_accession="GCF_000002285.1",
         circ_seqs=character(0)),

    list(assembly="CanFam2.0",
         date="2005/07/12",
         extra_info=c(breed="boxer"),
         assembly_accession="GCF_000002285.2",  # canFam2
         circ_seqs="MT"),

    list(assembly="CanFam3.1",
         date="2011/11/02",
         extra_info=c(breed="boxer"),
         assembly_accession="GCF_000002285.3",  # canFam3
         circ_seqs="MT"),

    list(assembly="UMICH_Zoey_3.1",
         date="2019/05/30",
         extra_info=c(breed="Great Dane"),
         assembly_accession="GCA_005444595.1",  # canFam5
         circ_seqs="chrM"),

    list(assembly="UU_Cfam_GSD_1.0",
         date="2020/03/10",
         extra_info=c(breed="German Shepherd"),
         assembly_accession="GCA_011100685.1",  # canFam4
         circ_seqs="chrM"),

    ## 147 sequences. Note that in December 2022, someone added a duplicated
    ## MT entry to GCF_000002285.5_Dog10K_Boxer_Tasha_assembly_report.txt,
    ## bringing the number of entries to 148, but without assigning a new
    ## RefSeq accession (still GCF_000002285.5), and without updating the
    ## assembly date either (still 2020/10/06). How can NCBI even allow this?
    ## Shocking!
    ## Anyways, on Jan 2, 2023, in GenomeInfoDb 1.34.6 (release) and 1.35.10
    ## (devel), I modified getChromInfoFromNCBI() so that it removes the bogus
    ## MT entry (which turns out to be the original one), bringing the number
    ## of entries back to 147. See .get_NCBI_chrom_info_from_accession() in
    ## R/getChromInfoFromNCBI.R for the details.
    list(assembly="Dog10K_Boxer_Tasha",
         date="2020/10/06",
         extra_info=c(breed="boxer"),
         assembly_accession="GCF_000002285.5",  # canFam6
         circ_seqs="MT")
)

