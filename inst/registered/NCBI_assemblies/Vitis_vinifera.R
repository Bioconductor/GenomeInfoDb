ORGANISM <- "Vitis vinifera"

### List of assemblies first by cultivar then by date.
### Yep, different genome assemblies can have the same name! (don't ask me why)
### Lookup by genome name will pick-up the first in the list.
ASSEMBLIES <- list(
    ## cultivar: Chardonnay
    list(assembly="AWRI_Vv-CHD_1.0",
         date="2019/01/10",
         extra_info=c(cultivar="Chardonnay"),
         assembly_accession="GCA_004011995.1",
         circ_seqs=character(0)),

    ## cultivar: PN40024
    list(assembly="8x_WGS",
         date="2007/09/19",
         extra_info=c(cultivar="PN40024"),
         assembly_accession="GCF_000003745.1",
         circ_seqs=character(0)),

    list(assembly="12X",
         date="2009/12/07",
         extra_info=c(cultivar="PN40024"),
         assembly_accession="GCF_000003745.2",
         circ_seqs=c("MT", "Pltd")),

    list(assembly="12X",
         date="2009/12/07",
         extra_info=c(cultivar="PN40024"),
         assembly_accession="GCF_000003745.3",
         circ_seqs=c("MT", "Pltd"))
)

