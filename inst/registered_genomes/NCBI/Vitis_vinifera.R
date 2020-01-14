ORGANISM <- "Vitis vinifera"

### List of assemblies first by cultivar then by date.
### Yep, different genome assemblies can have the same name! (don't ask me why)
### Lookup by genome name will pick-up the first in the list.
ASSEMBLIES <- list(
    ## cultivar: Chardonnay
    list(genome="AWRI_Vv-CHD_1.0",
         assembly_accession="GCA_004011995.1",
         infraspecific_name=c(cultivar="Chardonnay"),
         date="2019/01/10",
         circ_seqs=character(0)),

    ## cultivar: PN40024
    list(genome="8x_WGS",
         assembly_accession="GCF_000003745.1",
         infraspecific_name=c(cultivar="PN40024"),
         date="2007/09/19",
         circ_seqs=character(0)),

    list(genome="12X",
         assembly_accession="GCF_000003745.2",
         infraspecific_name=c(cultivar="PN40024"),
         date="2009/12/07",
         circ_seqs=c("MT", "Pltd")),

    list(genome="12X",
         assembly_accession="GCF_000003745.3",
         infraspecific_name=c(cultivar="PN40024"),
         date="2009/12/07",
         circ_seqs=c("MT", "Pltd"))
)

