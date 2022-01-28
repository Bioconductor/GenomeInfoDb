ORGANISM <- "Gallus gallus"

### List of assemblies by date.
### Yep, different genome assemblies can have the same name! (don't ask me why)
### Lookup by genome name will pick-up the first in the list.
ASSEMBLIES <- list(
    ## 17152 sequences.
    list(assembly="Gallus_gallus-2.1",
         date="2006/11/01",
         assembly_accession="GCF_000002315.2",  # galGal3
         circ_seqs="MT"),

    ## 17151 sequences, same as in the GCF_000002315.2 assembly above except
    ## that MT is missing.
    ## Also note that the 17151 sequences in GCF_000002315.1 are the same as
    ## in GCA_000002315.1 but 14 of them have different GenBank accessions.
    list(assembly="Gallus_gallus-2.1",
         date="2006/11/01",
         assembly_accession="GCF_000002315.1",
         circ_seqs=character(0)),

    list(assembly="Gallus_gallus-4.0",
         date="2011/11/22",
         assembly_accession="GCA_000002315.2",  # galGal4
         circ_seqs="MT"),

    list(assembly="Gallus_gallus-5.0",
         date="2015/12/16",
         assembly_accession="GCF_000002315.4",  # galGal5
         circ_seqs="MT"),

    list(assembly="Ogye1.0",
         date="2017/12/05",
         assembly_accession="GCA_002798355.1",
         circ_seqs="Ogye_chrM"),

    list(assembly="GRCg6",
         date="2018/02/02",
         assembly_accession="GCA_000002315.4",
         circ_seqs=character(0)),

    ## 464 sequences, same as in GCF_000002315.5.
    ## MT is NC_001323.1 (16775 bp), but this RefSeq record got removed
    ## by RefSeq stack, don't know why!
    list(assembly="GRCg6a",
         date="2018/03/27",
         assembly_accession="GCA_000002315.5",  # galGal6
         circ_seqs="MT"),

    ## 464 sequences, same as in the GCA_000002315.5/GCF_000002315.5 assembly
    ## above except that MT now is NC_040902.1 (16784 bp).
    list(assembly="GRCg6a",
         date="2018/03/27",
         assembly_accession="GCF_000002315.6",
         circ_seqs="MT"),

    ## 214 sequences, same as in GCA_016699485.1.
    list(assembly="bGalGal1.mat.broiler.GRCg7b",
         date="2021/01/19",
         assembly_accession="GCF_016699485.2",
         circ_seqs="MT"),

    ## 276 sequences.
    list(assembly="bGalGal1.pat.whiteleghornlayer.GRCg7w",
         date="2021/10/01",
         assembly_accession="GCA_016700215.2",
         circ_seqs=character(0)),

    ## 255 sequences, same as in GCA_016700215.1.
    ## Make sure to keep **after** the GCA_016700215.2 assembly above which
    ## has the same name but is more recent.
    list(assembly="bGalGal1.pat.whiteleghornlayer.GRCg7w",
         date="2021/01/19",
         assembly_accession="GCF_016700215.1",
         circ_seqs=character(0))
)

