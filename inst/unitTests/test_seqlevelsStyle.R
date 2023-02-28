test_basic <- function()
{
    filename <- system.file(package="GenomeInfoDb",  "extdata",
                            "dataFiles", "Homo_sapiens.txt")
    checkTrue(file.exists(filename))

    # check the format of the file
    data <- read.table(filename,header=TRUE,sep="\t")
    checkIdentical(c(25L, 7L), dim(data))
    checkIdentical(c('circular', 'auto', 'sex',
                     'NCBI', 'UCSC', 'dbSNP', 'Ensembl'),
                   colnames(data))

    #check if first 3 columns contain only true or false entries
    checkEquals(c(FALSE,TRUE),unique(data[,1]))
    checkEquals(c(TRUE,FALSE),unique(data[,2]))
    checkEquals(c(FALSE,TRUE),unique(data[,3]))
}

test_guessSpeciesStyle <- function()
{
    got  <- GenomeInfoDb:::.guessSpeciesStyle(c(paste0("chr",1:10)))
    checkEquals(unique(got$style), "UCSC")

    got <- GenomeInfoDb:::.guessSpeciesStyle(c(paste0("chr",1:22)))
    checkEquals(unique(got$style), "UCSC")

    got <- GenomeInfoDb:::.guessSpeciesStyle("chr2")
    checkEquals(unique(got$style), "UCSC")

    got <- GenomeInfoDb:::.guessSpeciesStyle("2")
    checkEquals(unique(got$style), c("NCBI","Ensembl","MSU6","JGI2.F","AGPvF"))

    got <- GenomeInfoDb:::.guessSpeciesStyle('T')
    checkEquals(unique(got$style), "JGI2.F")

    got <- GenomeInfoDb:::.guessSpeciesStyle(c("chr1","chr2","chr3",
        "chr1_gl000191_random", "chr1_gl000192_random"))
    checkEquals(unique(got$style), "UCSC")

    got <-  GenomeInfoDb:::.guessSpeciesStyle("h")
    checkEquals(got, NA)
}

test_seqlevelsStyle_character <- function()
{
    #1. correct seqnames
    got <- seqlevelsStyle(c(paste0('chr',1:20)))
    checkEquals(got,"UCSC")

    #2. mix seqnames from 2 styles for same organism
    got2 <- seqlevelsStyle(c('1','MT','Pltd','chr1'))
    checkEquals(got2,"NCBI")

    #3. mix seqnames from 2 different organisms
    got2 <- seqlevelsStyle(c('1','chr2RHet','chr3LHet'))
    checkEquals(got2,"UCSC")

    #4. incorrect seqnames
    checkException(seqlevelsStyle(c('234','567','acv')))

    #5. empty Seqinfo obj - with no seqnames
    checkException(seqlevelsStyle(c('')))
    checkException(seqlevelsStyle(GRanges()))

}

test_seqlevelsStyle_Seqinfo <- function()
{
    check_UCSC_NCBI_switch <- function(UCSC_genome, NCBI_assembly,
                                       nmapped, UCSC_nunmapped, NCBI_nunmapped)
    {
        ## Start with a Seqinfo object made from a UCSC genome.

        si1 <- Seqinfo(genome=UCSC_genome)
        UCSC_nseqlevels <- nmapped + UCSC_nunmapped
        checkEquals(UCSC_nseqlevels, length(si1))
        checkIdentical("UCSC", seqlevelsStyle(si1))

        si2 <- si1
        seqlevelsStyle(si2) <- "NCBI"
        ugenomes <- unique(genome(si2))
        ## UGLY HACK! We need to special-case hg38 because it contains 2
        ## sequences that do NOT belong to GRCh38.p14. But they can be
        ## found in GRCh38.p13!
        if (UCSC_genome == "hg38") {
            checkIdentical(c("GRCh38.p14", "GRCh38.p13"), ugenomes)
            checkIdentical("NCBI", seqlevelsStyle(si2))
	} else if (UCSC_nunmapped == 0L) {
            checkIdentical(NCBI_assembly, ugenomes)
            checkIdentical("NCBI", seqlevelsStyle(si2))
        } else {
            checkIdentical(c(NCBI_assembly, UCSC_genome), ugenomes)
            checkEquals(nmapped, sum(genome(si2) == NCBI_assembly))
            checkIdentical(c("NCBI", "UCSC"), seqlevelsStyle(si2))
        }

        seqlevelsStyle(si2) <- "UCSC"
        checkIdentical(si1, si2)

        ## Start with a Seqinfo object made from an NCBI assembly.

        si1 <- Seqinfo(genome=NCBI_assembly)
        NCBI_nseqlevels <- nmapped + NCBI_nunmapped
        checkEquals(NCBI_nseqlevels, length(si1))
        checkIdentical("NCBI", seqlevelsStyle(si1))

        si2 <- si1
        seqlevelsStyle(si2) <- "UCSC"
        ugenomes <- unique(genome(si2))
        if (NCBI_nunmapped <= 0L) {  # <= 0 because of hg38 special case!
            checkIdentical(UCSC_genome, ugenomes)
            checkIdentical("UCSC", seqlevelsStyle(si2))
        } else {
            ## 'ugenomes' will almost always be 'c(UCSC_genome, NCBI_assembly)'
            ## but the order is not 100% guaranteed (e.g. for
            ## musFur1/MusPutFur1.0 it's the opposite order).
            #checkIdentical(c(UCSC_genome, NCBI_assembly), ugenomes)
            checkEquals(2L, length(ugenomes))
            checkTrue(setequal(c(UCSC_genome, NCBI_assembly), ugenomes))
            checkEquals(nmapped, sum(genome(si2) == UCSC_genome))
            ## 'seqlevelsStyle(si2)' will almost always return c("UCSC", "NCBI")
            ## but the order is not 100% guaranteed (e.g. for
            ## musFur1/MusPutFur1.0 it's the opposite order).
            #checkIdentical(c("UCSC", "NCBI"), seqlevelsStyle(si2))
            checkEquals(2L, length(seqlevelsStyle(si2)))
            checkTrue(setequal(c("UCSC", "NCBI"), seqlevelsStyle(si2)))
        }

        seqlevelsStyle(si2) <- "NCBI"
        checkIdentical(si1, si2)
    }

    UCSC_NCBI <- list(
        # Field 1: UCSC_genome
        # Field 2: NCBI_assembly
        # Field 3: nb of seqlevels that are mapped between UCSC and NCBI
        # Field 4: nb of UCSC seqlevels that are not mapped to NCBI
        # Field 5: nb of NCBI seqlevels that are not mapped to UCSC
        #     1           2                               3      4        5
        list("apiMel2",  "Amel_2.0",                     16L,    1L,   7151L),
        #list("wuhCor1",  "ASM985889v3",                   1L,    0L,      0L),
        #list("bosTau6",  "Bos_taurus_UMD_3.1",         3317L,    0L,      0L),
        list("bosTau7",  "Btau_4.6.1",                11691L,    1L,      1L),
        list("bosTau8",  "Bos_taurus_UMD_3.1.1",       3179L,    0L,      0L),
        list("bosTau9",  "ARS-UCD1.2",                 2211L,    0L,      1L),
        list("calJac3",  "Callithrix jacchus-3.2",    14205L,    0L,      0L),
        list("calJac4",  "Callithrix_jacchus_cj1700_1.1", 964L,  0L,      0L),
        list("canFam3",  "CanFam3.1",                  3268L,    0L,      0L),
        list("canFam4",  "UU_Cfam_GSD_1.0",            2198L,    0L,      0L),
        list("canFam5",  "UMICH_Zoey_3.1",              794L,    0L,      0L),
        list("canFam6",  "Dog10K_Boxer_Tasha",          147L,    0L,      0L),
        list("ce6",      "WS190",                         7L,    0L,      0L),
        list("ce10",     "WBcel215",                      7L,    0L,      0L),
        list("ce11",     "WBcel235",                      7L,    0L,      0L),
        list("danRer7",  "Zv9",                        1133L,    0L,      0L),
        list("danRer10", "GRCz10",                     1061L,    0L,      0L),
        list("danRer11", "GRCz11",                     1923L,    0L,      0L),
        list("dm3",      "Release 5",                    14L,    1L,      0L),
        list("dm6",      "Release 6 plus ISO1 MT",     1870L,    0L,      0L),
        list("equCab2",  "EquCab2.0",                    33L,    1L,   9604L),
        list("equCab3",  "EquCab3.0",                  4701L,    0L,      0L),
        list("felCat9",  "Felis_catus_9.0",            4508L,    0L,      0L),
        list("galGal3",  "Gallus_gallus-2.1",            34L,   23L,  17118L),
        list("galGal4",  "Gallus_gallus-4.0",         15932L,    0L,      0L),
        list("galGal5",  "Gallus_gallus-5.0",         23475L,    0L,      0L),
        list("galGal6",  "GRCg6a",                      464L,    0L,      0L),
        list("hg15",     "NCBI33",                       24L,   20L,    140L),
        #list("hg16",     "NCBI34",                       24L,   18L,    138L),
        #list("hg17",     "NCBI35",                       26L,   20L,     86L),
        #list("hg18",     "NCBI36",                       26L,   23L,     97L),
        list("hg19",     "GRCh37.p13",                  297L,    1L,      0L),
        list("hg38",     "GRCh38.p14",                  711L,    0L,     -2L),
        list("hs1",      "T2T-CHM13v2.0",                25L,    0L,      0L),
        list("macFas5",  "Macaca_fascicularis_5.0",    7601L,    0L,      0L),
        list("mm8",      "MGSCv36",                      21L,   13L,    360L),
        list("mm9",      "MGSCv37",                      22L,   13L,    283L),
        list("mm10",     "GRCm38.p6",                   239L,    0L,      0L),
        list("mm39",     "GRCm39",                       61L,    0L,      0L),
        list("monDom5",  "MonDom5",                      10L,    1L,   9028L),
        list("musFur1",  "MusPutFur1.0",               7741L,    0L,     42L),
        list("panPan1",  "panpan1",                   10867L,    0L,      0L),
        list("panPan2",  "panpan1.1",                 10274L,    0L,      0L),
        list("panPan3",  "Mhudiblu_PPA_v0",            4293L,    0L,      0L),
        list("panTro2",  "Pan_troglodytes-2.1",          26L,   26L,  29214L),
        list("panTro3",  "Pan_troglodytes-2.1.3",     24131L,    1L,      0L),
        list("panTro4",  "Pan_troglodytes-2.1.4",     24129L,    0L,      0L),
        list("panTro5",  "Pan_tro 3.0",               44449L,    0L,      0L),
        list("panTro6",  "Clint_PTRv2",                4346L,    0L,      0L),
        list("rheMac2",  "Mmul_051212",                  21L,    1L, 122143L),
        list("rheMac3",  "CR_1.0",                    34102L,    1L,      0L),
        list("rheMac8",  "Mmul_8.0.1",               284728L,    0L,      0L),
        list("rheMac10", "Mmul_10",                    2939L,    0L,      0L),
        list("rn5",      "Rnor_5.0",                   2739L,    0L,      0L),
        list("rn6",      "Rnor_6.0",                    953L,    0L,      2L),
        list("rn7",      "mRatBN7.2",                   176L,    0L,      0L),
        list("sacCer3",  "R64",                          17L,    0L,      0L),
        list("susScr2",  "Sscrofa9.2",                   19L,    1L,      0L),
        list("susScr3",  "Sscrofa10.2",                4583L,    0L,      0L),
        list("susScr11", "Sscrofa11.1",                 613L,    0L,      0L),
        list("taeGut2",  "Taeniopygia_guttata-3.2.4", 37096L,    0L,      0L),
        list("xenLae2",  "Xenopus_laevis_v2",        108033L,    0L,      0L),
        list("xenTro10", "UCB_Xtro_10.0",               167L,    0L,      0L)
    )
    for (i in seq_along(UCSC_NCBI)) {
        args <- UCSC_NCBI[[i]]
        do.call(check_UCSC_NCBI_switch, args)
    }

    check_RefSeq_switch <- function(UCSC_genome, NCBI_assembly, UCSC_nunmapped)
    {
        is_RefSeq_accession <- GenomeInfoDb:::.is_RefSeq_accession

        ## UGLY HACK! Hardcode list of NCBI seqlevels (+ their UCSC names)
        ## NOT associated with a RefSeq accession in the "Full sequence
        ## report" for GRCh38.p14:
        if (NCBI_assembly == "GRCh38.p14") {
            problematic_NCBI_seqlevels <- c(
                "HSCHR10_1_CTG4",
                "HSCHR11_CTG1_UNLOCALIZED",
                "HSCHR22_UNLOCALIZED_CTG4",
                "HSCHRUN_RANDOM_CTG29"
            )
            ## Their UCSC names.
            problematic_UCSC_seqlevels <- c(
                "chr10_KI270825v1_alt",
                "chr11_KI270721v1_random",
                "chr22_KI270734v1_random",
                "chrUn_KI270752v1"
            )
        } else {
            problematic_NCBI_seqlevels <- character(0)
            problematic_UCSC_seqlevels <- character(0)
        }

        ## Start with a Seqinfo object made from a UCSC genome.

        si1 <- Seqinfo(genome=UCSC_genome)
        ## Remove problematic UCSC seqlevels:
        si1 <- si1[setdiff(seqlevels(si1), problematic_UCSC_seqlevels)]

        si2 <- si1
        seqlevelsStyle(si2) <- "NCBI"

        si3 <- si2
        seqlevelsStyle(si3) <- "RefSeq"
        checkIdentical(unname(genome(si2)), unname(genome(si3)))
        style <- seqlevelsStyle(si3)
        if (UCSC_nunmapped == 0L) {
            checkTrue(identical("RefSeq", style) ||
                      identical(c("RefSeq", "NCBI"), style))
        } else {
            checkIdentical(c("RefSeq", "UCSC"), style)
        }
        has_changed <- seqnames(si3) != seqnames(si2)
        ## UGLY HACK! We need to special-case hg38 because it contains 2
        ## sequences that do NOT belong to GRCh38.p14. But they can be
        ## found in GRCh38.p13!
        if (UCSC_genome == "hg38") {
            genome_is_ok <- genome(si2) %in% c("GRCh38.p14", "GRCh38.p13")
        } else {
            genome_is_ok <- genome(si2) == NCBI_assembly
        }
        checkTrue(all(genome_is_ok | !has_changed))
        checkTrue(all(is_RefSeq_accession(seqnames(si3)[has_changed])))

        si4 <- si1
        seqlevelsStyle(si4) <- "RefSeq"
        checkIdentical(si3, si4)
        seqlevelsStyle(si4) <- "UCSC"
        checkIdentical(si1, si4)

        seqlevelsStyle(si3) <- "NCBI"
        checkIdentical(si2, si3)

        seqlevelsStyle(si2) <- "UCSC"
        checkIdentical(si1, si2)

        ## Start with a Seqinfo object made from an NCBI assembly.

        si1 <- Seqinfo(genome=NCBI_assembly)
        ## Remove problematic NCBI seqlevels:
        si1 <- si1[setdiff(seqlevels(si1), problematic_NCBI_seqlevels)]

        si2 <- si1
        seqlevelsStyle(si2) <- "UCSC"

        si3 <- si1
        seqlevelsStyle(si3) <- "RefSeq"
        checkIdentical(unname(genome(si1)), unname(genome(si3)))
        style <- seqlevelsStyle(si3)
        checkTrue(identical("RefSeq", style) ||
                  identical(c("RefSeq", "NCBI"), style))
        has_changed <- seqnames(si3) != seqnames(si1)
        checkTrue(all(is_RefSeq_accession(seqnames(si3)[has_changed])))

        si4 <- si2
        seqlevelsStyle(si4) <- "RefSeq"
        checkIdentical(si3, si4)
        seqlevelsStyle(si4) <- "NCBI"
        checkIdentical(si1, si4)

        seqlevelsStyle(si2) <- "NCBI"
        checkIdentical(si1, si2)
    }

    ## Exclude some problematic genomes from the RefSeq switch check:
    ## - In the case of canFam5 and hs1, it's because chrM is not mapped
    ##   to a RefSeq accession, that is, RefSeqAccn is NA for chrM in
    ##   getChromInfoFromNCBI("UMICH_Zoey_3.1"), and for MT in
    ##   getChromInfoFromNCBI("T2T-CHM13v2.0").
    ## - For the other genomes (canFam4, rheMac3, panTro3), the reasons
    ##   causing check_RefSeq_switch() still need to be investigated!
    skip_RefSeq_switch <- c("canFam4", "canFam5", "hs1", "rheMac3", "panTro3")
    for (i in seq_along(UCSC_NCBI)) {
        args <- UCSC_NCBI[[i]][c(1L, 2L, 4L)]
        if (args[[1L]] %in% skip_RefSeq_switch)
            next
        do.call(check_RefSeq_switch, args)
    }
}

test_genomeStyles <- function()
{
    checkIdentical("data.frame", class(genomeStyles("Homo sapiens")))
    checkIdentical(c(25L, 7L), dim(genomeStyles("Homo sapiens")))
    checkException(genomeStyles("SAD"))
}

test_extractSeqlevels <- function()
{
    got <- extractSeqlevels("Homo sapiens", "UCSC" )
    checkEquals(25,length(got))
    checkEquals("character",class(got))

    checkException(extractSeqlevels("aaa","Homo sapiens"))
    checkException(extractSeqlevels("Drosophila melanogaster"))
}

test_extractSeqlevelsByGroup <- function()
{
    got <- extractSeqlevelsByGroup("Drosophila melanogaster","Ensembl","auto")
    checkEquals(5,length(got))
    checkEquals("character",class(got))

    checkException(extractSeqlevelsByGroup("aaa","Homo sapiens"))
    checkException(extractSeqlevelsByGroup("Drosophila melanogaster"))
    checkException(extractSeqlevelsByGroup("Homo sapiens","auto","NCBI"))
}

test_seqlevelsInGroup <- function()
{
    newch <- paste0("chr",c(1:22,"X","Y","M","1_gl000192_random","4_ctg9_hap1"))
    got1 <- seqlevelsInGroup(newch, group="sex")
    checkEquals(c("chrX","chrY"),got1)

    newchr <- as.character(c(1:22,"X","Y","MT"))
    got2 <- seqlevelsInGroup(newchr, group="all","Homo sapiens","NCBI")
    checkEquals(25,length(got2))
}

