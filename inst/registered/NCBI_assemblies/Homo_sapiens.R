ORGANISM <- "Homo sapiens"

### List of assemblies first by BioProject then by date.
ASSEMBLIES <- list(

    ## --- Human Genome Project ---
    ## BioProject: PRJNA31257

    list(assembly="NCBI33",
         assembly_level="Chromosome",
         date="2003/04/12",
         assembly_accession="GCF_000001405.8",   # hg15
         circ_seqs=character(0)),

    list(assembly="NCBI34",
         assembly_level="Chromosome",
         date="2004/02/04",
         assembly_accession="GCF_000001405.10",  # hg16
         circ_seqs=character(0)),

    list(assembly="NCBI35",
         assembly_level="Chromosome",
         date="2004/08/24",
         assembly_accession="GCF_000001405.11",  # hg17
         circ_seqs=character(0)),

    list(assembly="NCBI36",
         assembly_level="Chromosome",
         date="2006/03/03",
         assembly_accession="GCF_000001405.12",  # hg18
         circ_seqs=character(0)),

    list(assembly="GRCh37",
         assembly_level="Chromosome",
         date="2009/02/27",
         assembly_accession="GCF_000001405.13",
         circ_seqs=character(0)),

    list(assembly="GRCh37.p13",
         assembly_level="Chromosome",
         date="2013/06/28",
         assembly_accession="GCF_000001405.25",  # hg19
         circ_seqs="MT"),

    ## 455 sequences.
    list(assembly="GRCh38",
         assembly_level="Chromosome",
         date="2013/12/17",
         assembly_accession="GCF_000001405.26",
         circ_seqs="MT"),

    ## 471 sequences.
    list(assembly="GRCh38.p1",
         assembly_level="Chromosome",
         date="2014/10/03",
         assembly_accession="GCF_000001405.27",
         circ_seqs="MT"),

    ## 486 sequences.
    list(assembly="GRCh38.p2",
         assembly_level="Chromosome",
         date="2014/12/05",
         assembly_accession="GCF_000001405.28",
         circ_seqs="MT"),

    ## 494 sequences.
    list(assembly="GRCh38.p3",
         assembly_level="Chromosome",
         date="2015/04/03",
         assembly_accession="GCF_000001405.29",
         circ_seqs="MT"),

    ## 510 sequences.
    list(assembly="GRCh38.p4",
         assembly_level="Chromosome",
         date="2015/06/25",
         assembly_accession="GCF_000001405.30",
         circ_seqs="MT"),

    ## 517 sequences.
    list(assembly="GRCh38.p5",
         assembly_level="Chromosome",
         date="2015/09/22",
         assembly_accession="GCF_000001405.31",
         circ_seqs="MT"),

    ## 521 sequences.
    list(assembly="GRCh38.p6",
         assembly_level="Chromosome",
         date="2015/12/21",
         assembly_accession="GCF_000001405.32",
         circ_seqs="MT"),

    ## 525 sequences.
    list(assembly="GRCh38.p7",
         assembly_level="Chromosome",
         date="2016/03/21",
         assembly_accession="GCF_000001405.33",
         circ_seqs="MT"),

    ## 543 sequences.
    list(assembly="GRCh38.p8",
         assembly_level="Chromosome",
         date="2016/06/30",
         assembly_accession="GCF_000001405.34",
         circ_seqs="MT"),

    ## 551 sequences.
    list(assembly="GRCh38.p9",
         assembly_level="Chromosome",
         date="2016/09/26",
         assembly_accession="GCF_000001405.35",
         circ_seqs="MT"),

    ## 557 sequences.
    list(assembly="GRCh38.p10",
         assembly_level="Chromosome",
         date="2017/01/06",
         assembly_accession="GCF_000001405.36",
         circ_seqs="MT"),

    ## 578 sequences.
    list(assembly="GRCh38.p11",
         assembly_level="Chromosome",
         date="2017/06/14",
         assembly_accession="GCF_000001405.37",
         circ_seqs="MT"),

    ## 595 sequences.
    list(assembly="GRCh38.p12",
         assembly_level="Chromosome",
         date="2017/12/21",
         assembly_accession="GCF_000001405.38",
         circ_seqs="MT"),

    ## 640 sequences.
    list(assembly="GRCh38.p13",
         assembly_level="Chromosome",
         date="2019/02/28",
         ## Ensembl uses GCA_000001405.28 for homo_sapiens in release
         ## 99 which is the same as GCF_000001405.39 **except** that
         ## the former does NOT have the RefSeq accessions!
         assembly_accession="GCF_000001405.39",  # hg38
         circ_seqs="MT"),

    ## 709 sequences.
    list(assembly="GRCh38.p14",
         assembly_level="Chromosome",
         date="2022/02/03",
         assembly_accession="GCF_000001405.40",
         circ_seqs="MT"),

    ## --- T2T-CHM13 project by the Telomere-to-Telomere Consortium ---
    ## BioProject: PRJNA559484

    ## 359 sequences.
    list(assembly="T2T-CHM13v0.7",
         assembly_level="Chromosome",
         date="2020/01/22",
         assembly_accession="GCA_009914755.1",
         circ_seqs=character(0)),

    ## 24 sequences (Y missing).
    list(assembly="T2T-CHM13v1.0",
         assembly_level="Chromosome",
         date="2021/01/26",
         assembly_accession="GCA_009914755.2",
         circ_seqs="MT"),

    ## 24 sequences (Y missing).
    list(assembly="T2T-CHM13v1.1",
         assembly_level="Complete Genome",
         date="2021/05/07",
         assembly_accession="GCA_009914755.3",
         circ_seqs="MT"),

    ## 25 sequences.
    list(assembly="T2T-CHM13v2.0",
         assembly_level="Complete Genome",
         date="2022/01/24",
         assembly_accession="GCA_009914755.4",   # hs1
         circ_seqs="MT")
)

