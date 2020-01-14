### Should work as a standlone, self-contained script.
### Must define at least:
###   o GENOME:              Single non-empty string.
###   o ORGANISM:            Single non-empty string.
###   o ASSEMBLED_MOLECULES: Character vector with no NAs, no empty strings,
###                          and no duplicates.
### Can also define:
###   o CIRC_SEQS:           Character vector (subset of ASSEMBLED_MOLECULES).
###   o GET_CHROM_SIZES:     Function with 1 argument. Must return a 2-column
###                          data.frame with columns "chrom" and "size".
###   o NCBI_LINKER:         Named list.
GENOME <- "hg19"
ORGANISM <- "Homo sapiens"
ASSEMBLED_MOLECULES <- paste0("chr", c(1:22, "X", "Y", "M"))

CIRC_SEQS <- "chrM"

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
    m33 <- m34 <- match(m3[ , 3L], c(paste0("hap", 1:7), "random"))
    stopifnot(!anyNA(m33))
    m33[m33 <= 7L] <- 1L
    oo3 <- order(m33, m31, m34, m3[ , 2L])
    idx3 <- idx3[oo3]

    idx2 <- which(npart == 2L)
    m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
    stopifnot(all(m2[ , 1L] == "chrUn"))
    oo2 <- order(m2[ , 2L])
    idx2 <- idx2[oo2]

    c(idx1, idx3, idx2)
}

GET_CHROM_SIZES <- function(goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    chrom_sizes <- GenomeInfoDb:::fetch_chrom_sizes_from_UCSC(GENOME,
                                              goldenPath.url=goldenPath.url)
    oo <- .order_seqlevels(chrom_sizes[ , "chrom"])
    S4Vectors:::extract_data_frame_rows(chrom_sizes, oo)
}

### Valid NCBI_LINKER components:
### - assembly_accession: single non-empty string.
### - AssemblyUnits: character vector.
### - special_mappings: named character vector.
### - unmapped_seqs: named list of character vectors.
### - drop_unmapped: TRUE or FALSE.
NCBI_LINKER <- list(
    assembly_accession="GCF_000001405.13",
    ## Special renaming of the 9 alternate scaffolds:
    special_mappings=c(chr4_ctg9_hap1="HSCHR4_1_CTG9",
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

