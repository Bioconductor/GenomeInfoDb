ORGANISM <- "Bos taurus"

### List of assemblies by date.
### Yep, different genome assemblies can have the same name!
### Lookup by genome name will pick-up the first in the list.
ASSEMBLIES <- list(
    list(genome="Bos_taurus_UMD_3.1",
         assembly_accession="GCF_000003055.4",  # bosTau6
         date="2009/12/01",
         circ_seqs="MT"),

    list(genome="Bos_taurus_UMD_3.1.1",
         assembly_accession="GCF_000003055.5",  # bosTau8
         date="2009/12/01",
         circ_seqs="MT"),

    list(genome="Btau_4.6.1",
         assembly_accession="GCA_000003205.4",  # bosTau7
         date="2011/11/02",
         circ_seqs=character(0)),

    list(genome="Bos_taurus_UMD_3.1.1",
         assembly_accession="GCA_000003055.5",
         date="2014/11/25",
         circ_seqs="MT"),

    list(genome="Btau_5.0.1",
         assembly_accession="GCA_000003205.6",
         date="2015/11/19",
         circ_seqs=character(0)),

    list(genome="ARS-UCD1.2",
         assembly_accession="GCA_002263795.2",  # bosTau9
         date="2018/04/11",
         circ_seqs=c("MT", "M"))
)

