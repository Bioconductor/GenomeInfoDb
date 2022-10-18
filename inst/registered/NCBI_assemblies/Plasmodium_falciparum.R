ORGANISM <- "Plasmodium falciparum"

### List of assemblies first by isolate, then by date.
ASSEMBLIES <- list(
    ## --- isolate: 3D7 (NCBI:txid36329) ---

    ## 15 sequences.
    list(assembly="ASM171558v1",
         assembly_level="Chromosome",
         date="2016/08/29",
         extra_info=c(isolate="3D7"),
         assembly_accession="GCA_001715585.1",
         ## This assembly includes 1 non-nlucear sequence:
         ##   - MT: Mitochondrion, GenBank-Accn CP017005.1
         ## According to https://www.ncbi.nlm.nih.gov/nuccore/CP017005.1,
         ## this sequence is circular.
         circ_seqs="MT"),

    ## 16 sequences.
    list(assembly="GCA_000002765",
         assembly_level="Complete Genome",
         date="2019/10/22",
         extra_info=c(isolate="3D7"),
         assembly_accession="GCA_000002765.3",
         ## This assembly includes 2 non-nlucear sequences:
         ##   - API: Apicoplast, GenBank-Accn LR605956.1
         ##   - MIT: Mitochondrion, GenBank-Accn LR605957.1
         ## According to https://www.ncbi.nlm.nih.gov/nuccore/LR605956.1 and
         ## https://www.ncbi.nlm.nih.gov/nuccore/LR605957.1, none of them is
         ## circular!
         circ_seqs=character(0)),

    ## --- isolate: HB3 (NCBI:txid137071) ---

    ## 1189 sequences (Primary Assembly unit only).
    list(assembly="ASM14966v1",
         assembly_level="Scaffold",
         date="2009/09/22",
         extra_info=c(isolate="HB3"),
         assembly_accession="GCA_000149665.1",
         circ_seqs=character(0)),

    ## 1191 sequences (same as ASM14966v1 + non-nuclear assembly unit).
    list(assembly="ASM14966v2",
         assembly_level="Scaffold",
         date="2009/09/22",
         extra_info=c(isolate="HB3"),
         assembly_accession="GCA_000149665.2",
         ## This assembly includes 2 non-nlucear sequences:
         ##   - MT: Mitochondrion, GenBank-Accn DQ642845.1
         ##   - Pltd: Apicoplast, GenBank-Accn DQ642846.1
         ## According to https://www.ncbi.nlm.nih.gov/nuccore/DQ642845.1 and
         ## https://www.ncbi.nlm.nih.gov/nuccore/DQ642846.1, both of them are
         ## circular.
         circ_seqs=c("MT", "Pltd")),

    ## 28 sequences.
    list(assembly="PfHB3-3",
         assembly_level="Chromosome",
         date="2018/12/11",
         extra_info=c(isolate="HB3"),
         assembly_accession="GCA_900631985.1",
         ## This assembly includes 2 non-nlucear sequences:
         ##   - PfHB3_MT: Mitochondrion, GenBank-Accn LR131352.1
         ##   - PfHB3_API: Apicoplast, GenBank-Accn LR131353.1
         ## According to https://www.ncbi.nlm.nih.gov/nuccore/LR131352.1 and
         ## https://www.ncbi.nlm.nih.gov/nuccore/LR131353.1, none of them is
         ## circular!
         circ_seqs=character(0))
)

