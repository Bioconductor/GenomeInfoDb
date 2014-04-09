### =========================================================================
### Some low-level (non exported) utility functions.
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetch_GRCh38_hg38_seqlevels()
###

fetch_GenBankAccn2seqlevel_for_GRCh38 <- function()
{
    NCBI_assembly_report_url <- paste(
        "ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes",
        "Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38",
        "GCA_000001405.15_GRCh38_assembly_report.txt", sep="/")
    colnames <- c("Sequence-Name", "Sequence-Role", "Assigned-Molecule",
                  "Assigned-Molecule-Location/Type", "GenBank-Accn",
                  "Relationship", "RefSeq-Accn", "Assembly-Unit")
    NCBI_assembly_report <- read.table(NCBI_assembly_report_url, sep="\t",
                                       col.names=colnames,
                                       stringsAsFactors=FALSE)
    GenBank_accn <- NCBI_assembly_report[[5L]]
    stopifnot(!any(duplicated(GenBank_accn)))
    ans <- NCBI_assembly_report[[1L]]
    names(ans) <- GenBank_accn
    ans
}

fetch_GenBankAccn2seqlevel_for_hg38 <- function(GRCh38_accn2seqlevel)
{
    suppressMessages(require(IRanges, quietly=TRUE))  # for elementLengths()
    hg38_chrominfo_url <- paste(
        "http://hgdownload.soe.ucsc.edu/goldenPath",
        "hg38/database/chromInfo.txt.gz", sep="/")
    destfile <- tempfile()
    download.file(hg38_chrominfo_url, destfile, quiet=TRUE)
    colnames <- c("chrom", "size", "fileName")
    chrominfo <- read.table(destfile, sep="\t", quote="",
                            col.names=colnames, comment.char="",
                            stringsAsFactors=FALSE)
    hg38_seqlevels <- chrominfo$chrom
    tmp <- lapply(names(GRCh38_accn2seqlevel), grep, hg38_seqlevels)
    idx1 <- which(elementLengths(tmp) == 1L)
    tmp[idx1] <- as.list(hg38_seqlevels[unlist(tmp[idx1], use.names=FALSE)])
    NCBI_main_chromosomes <- c(1:22, "X", "Y", "MT")
    UCSC_main_chromosomes <- paste0("chr", c(1:22, "X", "Y", "M"))
    idx0 <- which(elementLengths(tmp) == 0L)
    stopifnot(all(GRCh38_accn2seqlevel[idx0] == NCBI_main_chromosomes))
    tmp[idx0] <- as.list(UCSC_main_chromosomes)
    stopifnot(all(elementLengths(tmp) == 1L))
    ans <- unlist(tmp, use.names=FALSE)
    names(ans) <- names(GRCh38_accn2seqlevel)
    ans
}

fetch_GRCh38_hg38_seqlevels <- function()
{
    GRCh38_accn2seqlevel <- fetch_GenBankAccn2seqlevel_for_GRCh38()
    hg38_accn2seqlevel <-
        fetch_GenBankAccn2seqlevel_for_hg38(GRCh38_accn2seqlevel)
    data.frame(accn=names(GRCh38_accn2seqlevel),
               GRCh38_seqlevel=unname(GRCh38_accn2seqlevel),
               hg38_seqlevel=unname(hg38_accn2seqlevel),
               stringsAsFactors=FALSE)
}

