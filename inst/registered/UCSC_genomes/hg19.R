GENOME <- "hg19"
ORGANISM <- "Homo sapiens"
## MT only added to hg19 in March 2020!
## See https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome-announce/whsRXqNsn24
ASSEMBLED_MOLECULES <- paste0("chr", c(1:22, "X", "Y", "M", "MT"))
CIRC_SEQS <- c("chrM", "chrMT")

library(IRanges)       # for CharacterList()
library(GenomeInfoDb)  # for fetch_chrom_sizes_from_UCSC()

.order_seqlevels <- function(seqlevels)
{
    tmp <- CharacterList(strsplit(seqlevels, "_"))
    npart <- lengths(tmp)
    stopifnot(all(npart <= 3L))

    idx1 <- which(npart == 1L)
    stopifnot(length(idx1) == length(ASSEMBLED_MOLECULES))
    oo1 <- match(ASSEMBLED_MOLECULES, seqlevels[idx1])
    stopifnot(!anyNA(oo1))
    idx1 <- idx1[oo1]

    idx3 <- which(npart == 3L)
    m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
    m31 <- match(m3[ , 1L], ASSEMBLED_MOLECULES)
    stopifnot(!anyNA(m31))
    m33 <- match(m3[ , 3L], c(paste0("hap", 1:7), "random", "fix", "alt"))
    stopifnot(!anyNA(m33))
    m33b <- m33
    m33[m33 <= 7L] <- 1L
    oo3 <- order(m33, m31, m33b, m3[ , 2L])
    idx3 <- idx3[oo3]
    idx3B_len <- sum(m33 >= 9L)
    if (length(idx3B_len) == 0L) {
        idx3A <- idx3
        idx3B <- integer(0)
    } else {
        idx3A <- head(idx3, n=-idx3B_len)
        idx3B <- tail(idx3, n=idx3B_len)
    }

    idx2 <- which(npart == 2L)
    m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
    stopifnot(all(m2[ , 1L] == "chrUn"))
    oo2 <- order(m2[ , 2L])
    idx2 <- idx2[oo2]

    c(idx1, idx3A, idx2, idx3B)
}

FETCH_ORDERED_CHROM_SIZES <-
    function(goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    chrom_sizes <- GenomeInfoDb:::fetch_chrom_sizes_from_UCSC(GENOME,
                                              goldenPath.url=goldenPath.url)
    oo <- .order_seqlevels(chrom_sizes[ , "chrom"])
    S4Vectors:::extract_data_frame_rows(chrom_sizes, oo)
}

NCBI_LINKER <- list(
    ## Used to be mapped to GCF_000001405.13 (GRCh37) but is now mapped
    ## to GRCh37.p13 (since March 2020, see link to announcement on
    ## genome-announce mailing list at top of this file).
    assembly_accession="GCF_000001405.25",
    ## According to the UCSC-style-name column of the assembly report for
    ## GRCh37.p13, MT in GRCh37.p13 is mapped to chrM in hg19, which is wrong!
    ## The truth is that it's mapped to the new chrMT sequence in hg19 (chrMT
    ## was only added to hg19 in March 2020, see comment at top of this file).
    ## So we use the 'special_mappings' field to explicitly map chrMT to MT.
    special_mappings=c(chrMT="MT",
                       ## Special renaming of the 9 alternate scaffolds:
                       chr4_ctg9_hap1="HSCHR4_1_CTG9",
                       chr6_apd_hap1="HSCHR6_MHC_APD_CTG1",
                       chr6_cox_hap2="HSCHR6_MHC_COX_CTG1",
                       chr6_dbb_hap3="HSCHR6_MHC_DBB_CTG1",
                       chr6_mann_hap4="HSCHR6_MHC_MANN_CTG1",
                       chr6_mcf_hap5="HSCHR6_MHC_MCF_CTG1",
                       chr6_qbl_hap6="HSCHR6_MHC_QBL_CTG1",
                       chr6_ssto_hap7="HSCHR6_MHC_SSTO_CTG1",
                       chr17_ctg5_hap1="HSCHR17_1_CTG5"),
    unmapped_seqs=list(`assembled-molecule`="chrM")
)

ENSEMBL_LINKER <- "ucscToEnsembl"

