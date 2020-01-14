ORGANISM <- "Bos taurus"

### List of assemblies first by breed then by date.
### Yep, different genome assemblies can have the same name! (don't ask me why)
### Lookup by genome name will pick-up the first in the list.
ASSEMBLIES <- list(
    ## breed: Hereford
    list(genome="UMD Bos_taurus 2.0",
         assembly_accession="GCA_000003055.1",
         infraspecific_name=c(breed="Hereford"),
         date="2009/04/30",
         circ_seqs=character(0)),

    list(genome="Bos_taurus_UMD_3.0",
         assembly_accession="GCA_000003055.2",
         infraspecific_name=c(breed="Hereford"),
         date="2009/09/09",
         circ_seqs=character(0)),

    list(genome="Bos_taurus_UMD_3.1",
         assembly_accession="GCF_000003055.3",
         infraspecific_name=c(breed="Hereford"),
         date="2009/12/01",
         circ_seqs=character(0)),

    list(genome="Bos_taurus_UMD_3.1",
         assembly_accession="GCF_000003055.4",  # bosTau6
         infraspecific_name=c(breed="Hereford"),
         date="2009/12/01",
         circ_seqs="MT"),

    list(genome="Bos_taurus_UMD_3.1.1",
         assembly_accession="GCF_000003055.5",  # bosTau8
         infraspecific_name=c(breed="Hereford"),
         date="2009/12/01",
         circ_seqs="MT"),

    list(genome="Btau_4.2",
         assembly_accession="GCF_000003205.3",
         infraspecific_name=c(breed="Hereford"),
         date="2011/05/12",
         circ_seqs=character(0)),

    list(genome="Btau_4.2",
         assembly_accession="GCF_000003205.4",
         infraspecific_name=c(breed="Hereford"),
         date="2011/05/12",
         circ_seqs=character(0)),

    list(genome="Btau_4.6",
         assembly_accession="GCA_000003205.3",
         infraspecific_name=c(breed="Hereford"),
         date="2011/09/14",
         circ_seqs=character(0)),

    list(genome="Btau_4.6.1",
         assembly_accession="GCA_000003205.4",  # bosTau7
         infraspecific_name=c(breed="Hereford"),
         date="2011/11/02",
         circ_seqs=character(0)),

    list(genome="Bos_taurus_UMD_3.1.1",
         assembly_accession="GCA_000003055.5",
         infraspecific_name=c(breed="Hereford"),
         date="2014/11/25",
         circ_seqs="MT"),

    list(genome="Btau_5.0",
         assembly_accession="GCA_000003205.5",
         infraspecific_name=c(breed="Hereford"),
         date="2015/09/30",
         circ_seqs=character(0)),

    list(genome="Btau_5.0.1",
         assembly_accession="GCA_000003205.6",
         infraspecific_name=c(breed="Hereford"),
         date="2015/11/19",
         circ_seqs=character(0)),

    list(genome="ARS-UCD1.2",
         ## GCA_002263795.2 is the same as GCF_002263795.1 but the former
         ## has the UCSC style names.
         assembly_accession="GCA_002263795.2",  # bosTau9
         infraspecific_name=c(breed="Hereford"),
         date="2018/04/11",
         circ_seqs=c("MT", "M"))
)

