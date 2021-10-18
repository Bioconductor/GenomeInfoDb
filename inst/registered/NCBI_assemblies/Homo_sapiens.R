ORGANISM <- "Homo sapiens"

### List of assemblies by date.
ASSEMBLIES <- list(
    list(assembly="NCBI33",
         date="2003/04/12",
         assembly_accession="GCF_000001405.8",   # hg15
         circ_seqs=character(0)),

    list(assembly="NCBI34",
         date="2004/02/04",
         assembly_accession="GCF_000001405.10",  # hg16
         circ_seqs=character(0)),

    list(assembly="NCBI35",
         date="2004/08/24",
         assembly_accession="GCF_000001405.11",  # hg17
         circ_seqs=character(0)),

    list(assembly="NCBI36",
         date="2006/03/03",
         assembly_accession="GCF_000001405.12",  # hg18
         circ_seqs=character(0)),

    list(assembly="GRCh37",
         date="2009/02/27",
         assembly_accession="GCF_000001405.13",
         circ_seqs=character(0)),

    list(assembly="GRCh37.p13",
         date="2013/06/28",
         assembly_accession="GCF_000001405.25",  # hg19
         circ_seqs="MT"),

    list(assembly="GRCh38",
         date="2013/12/17",
         assembly_accession="GCF_000001405.26",
         circ_seqs="MT"),

    list(assembly="GRCh38.p1",
         date="2014/10/03",
         assembly_accession="GCF_000001405.27",
         circ_seqs="MT"),

    list(assembly="GRCh38.p2",
         date="2014/12/05",
         assembly_accession="GCF_000001405.28",
         circ_seqs="MT"),

    list(assembly="GRCh38.p3",
         date="2015/04/03",
         assembly_accession="GCF_000001405.29",
         circ_seqs="MT"),

    list(assembly="GRCh38.p4",
         date="2015/06/25",
         assembly_accession="GCF_000001405.30",
         circ_seqs="MT"),

    list(assembly="GRCh38.p5",
         date="2015/09/22",
         assembly_accession="GCF_000001405.31",
         circ_seqs="MT"),

    list(assembly="GRCh38.p6",
         date="2015/12/21",
         assembly_accession="GCF_000001405.32",
         circ_seqs="MT"),

    list(assembly="GRCh38.p7",
         date="2016/03/21",
         assembly_accession="GCF_000001405.33",
         circ_seqs="MT"),

    list(assembly="GRCh38.p8",
         date="2016/06/30",
         assembly_accession="GCF_000001405.34",
         circ_seqs="MT"),

    list(assembly="GRCh38.p9",
         date="2016/09/26",
         assembly_accession="GCF_000001405.35",
         circ_seqs="MT"),

    list(assembly="GRCh38.p10",
         date="2017/01/06",
         assembly_accession="GCF_000001405.36",
         circ_seqs="MT"),

    list(assembly="GRCh38.p11",
         date="2017/06/14",
         assembly_accession="GCF_000001405.37",
         circ_seqs="MT"),

    list(assembly="GRCh38.p12",
         date="2017/12/21",
         assembly_accession="GCF_000001405.38",
         circ_seqs="MT"),

    list(assembly="GRCh38.p13",
         date="2019/02/28",
         ## Ensembl uses GCA_000001405.28 for homo_sapiens in release
         ## 99 which is the same as GCF_000001405.39 **except** that
         ## the former does NOT have the RefSeq accessions!
         assembly_accession="GCF_000001405.39",  # hg38
         circ_seqs="MT")
)

