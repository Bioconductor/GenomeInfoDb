ORGANISM <- "Bos taurus"

### List of assemblies first by breed then by date.
### Yep, different genome assemblies can have the same name! (don't ask me why)
### Lookup by genome name will pick-up the first in the list.
ASSEMBLIES <- list(
    ## breed: Hereford
    list(assembly="UMD Bos_taurus 2.0",
         assembly_level="Chromosome",
         date="2009/04/30",
         extra_info=c(breed="Hereford"),
         assembly_accession="GCA_000003055.1",
         circ_seqs=character(0)),

    list(assembly="Bos_taurus_UMD_3.0",
         assembly_level="Chromosome",
         date="2009/09/09",
         extra_info=c(breed="Hereford"),
         assembly_accession="GCA_000003055.2",
         circ_seqs=character(0)),

    list(assembly="Bos_taurus_UMD_3.1",
         assembly_level="Chromosome",
         date="2009/12/01",
         extra_info=c(breed="Hereford"),
         assembly_accession="GCF_000003055.3",
         circ_seqs=character(0)),

    list(assembly="Bos_taurus_UMD_3.1",
         assembly_level="Chromosome",
         date="2009/12/01",
         extra_info=c(breed="Hereford"),
         assembly_accession="GCF_000003055.4",  # bosTau6
         circ_seqs="MT"),

    list(assembly="Bos_taurus_UMD_3.1.1",
         assembly_level="Chromosome",
         date="2009/12/01",
         extra_info=c(breed="Hereford"),
         assembly_accession="GCF_000003055.5",  # bosTau8
         circ_seqs="MT"),

    list(assembly="Btau_4.2",
         assembly_level="Chromosome",
         date="2011/05/12",
         extra_info=c(breed="Hereford"),
         assembly_accession="GCF_000003205.3",
         circ_seqs=character(0)),

    list(assembly="Btau_4.2",
         assembly_level="Chromosome",
         date="2011/05/12",
         extra_info=c(breed="Hereford"),
         assembly_accession="GCF_000003205.4",
         circ_seqs=character(0)),

    list(assembly="Btau_4.6",
         assembly_level="Chromosome",
         date="2011/09/14",
         extra_info=c(breed="Hereford"),
         assembly_accession="GCA_000003205.3",
         circ_seqs=character(0)),

    list(assembly="Btau_4.6.1",
         assembly_level="Chromosome",
         date="2011/11/02",
         extra_info=c(breed="Hereford"),
         assembly_accession="GCA_000003205.4",  # bosTau7
         circ_seqs=character(0)),

    list(assembly="Bos_taurus_UMD_3.1.1",
         assembly_level="Chromosome",
         date="2014/11/25",
         extra_info=c(breed="Hereford"),
         assembly_accession="GCA_000003055.5",
         circ_seqs="MT"),

    list(assembly="Btau_5.0",
         assembly_level="Chromosome",
         date="2015/09/30",
         extra_info=c(breed="Hereford"),
         assembly_accession="GCA_000003205.5",
         circ_seqs=character(0)),

    list(assembly="Btau_5.0.1",
         assembly_level="Chromosome",
         date="2015/11/19",
         extra_info=c(breed="Hereford"),
         assembly_accession="GCA_000003205.6",
         circ_seqs=character(0)),

    list(assembly="ARS-UCD1.2",
         assembly_level="Chromosome",
         date="2018/04/11",
         extra_info=c(breed="Hereford"),
         ## GCA_002263795.2 is the same as GCF_002263795.1 but the former
         ## has the UCSC style names.
         assembly_accession="GCA_002263795.2",  # bosTau9
         circ_seqs=c("MT", "M"))
)

