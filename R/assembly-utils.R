### =========================================================================
### Some low-level (non exported) utility functions.
### -------------------------------------------------------------------------


.NCBI_ASSEMBLY_REPORTS_URL <-
    "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/"

.GENBANK_ASSEMBLY_ACCESSION_PREFIX <- "GCA_"
.TOP_LEVEL_GCA_URL <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/"
### Paths: genomes/all/GCA/xxx/xxx/xxx/

.REFSEQ_ASSEMBLY_ACCESSION_PREFIX <- "GCF_"
.TOP_LEVEL_GCF_URL <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
### Paths: genomes/all/GCF/xxx/xxx/xxx/

.is_genbank_assembly_accession <- function(x)
{
    ## We use %in% instead of == to be NA-proof.
    substr(x, 1L, nchar(.GENBANK_ASSEMBLY_ACCESSION_PREFIX)) %in%
        .GENBANK_ASSEMBLY_ACCESSION_PREFIX
}

.is_refseq_assembly_accession <- function(x)
{
    ## We use %in% instead of == to be NA-proof.
    substr(x, 1L, nchar(.REFSEQ_ASSEMBLY_ACCESSION_PREFIX)) %in%
        .REFSEQ_ASSEMBLY_ACCESSION_PREFIX
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetch_assembly_summary()
###

### Performs some quick sanity checks on the assembly summary.
.check_assembly_summary <- function(assembly_summary, genbank_or_refseq)
{
    ## Check "assembly_accession" prefixes.
    assembly_accession <- assembly_summary[ , "assembly_accession"]
    if (genbank_or_refseq == "genbank") {
        is_valid_accession <- .is_genbank_assembly_accession
    } else {
        is_valid_accession <- .is_refseq_assembly_accession
    }
    stopifnot(all(is_valid_accession(assembly_accession)))

    ## Check "assembly_accession" uniqueness.
    stopifnot(anyDuplicated(assembly_accession) == 0L)

    ## Check "paired_asm_comp" levels.
    paired_asm_comp <- assembly_summary[ , "paired_asm_comp"]
    paired_asm_comp_levels <- c("identical", "different", NA)
    stopifnot(all(paired_asm_comp %in% paired_asm_comp_levels))

    ## Check gbrs_paired_asm/paired_asm_comp consistency.
    gbrs_paired_asm <- assembly_summary[ , "gbrs_paired_asm"]
    idx1 <- which(is.na(gbrs_paired_asm))
    idx2 <- which(is.na(paired_asm_comp))
    stopifnot(identical(idx1, idx2))
    if (length(idx1) != 0L)
        gbrs_paired_asm <- gbrs_paired_asm[-idx1]

    ## Check "gbrs_paired_asm" prefixes.
    if (genbank_or_refseq == "genbank") {
        is_valid_accession <- .is_refseq_assembly_accession
    } else {
        is_valid_accession <- .is_genbank_assembly_accession
    }
    stopifnot(all(is_valid_accession(gbrs_paired_asm)))

    ## Check "gbrs_paired_asm" uniqueness.
    stopifnot(anyDuplicated(gbrs_paired_asm) == 0L)
}

.assembly_summary_cache <- new.env(parent=emptyenv())

fetch_assembly_summary <- function(genbank_or_refseq, quiet=FALSE)
{
    objname <- paste0(genbank_or_refseq, "_assembly_summary")
    ans <- try(get(objname, envir=.assembly_summary_cache, inherits=FALSE),
               silent=TRUE)
    if (!is(ans, "try-error"))
        return(ans)
    assembly_summary_filename <- paste0(objname, ".txt")
    url <- paste0(.NCBI_ASSEMBLY_REPORTS_URL, assembly_summary_filename)
    destfile <- tempfile()
    download.file(url, destfile, quiet=quiet)
    ans <- read.table(destfile, header=TRUE, sep="\t", quote="",
                                na.strings="na", skip=1L, comment.char="",
                                stringsAsFactors=FALSE)
    expected_colnames <- c(
        "X..assembly_accession", "bioproject", "biosample", "wgs_master",
        "refseq_category", "taxid", "species_taxid", "organism_name",
        "infraspecific_name", "isolate", "version_status", "assembly_level",
        "release_type", "genome_rep", "seq_rel_date", "asm_name", "submitter",
        "gbrs_paired_asm", "paired_asm_comp", "ftp_path",
        "excluded_from_refseq")
    if (!identical(expected_colnames, colnames(ans)))
        stop(url, " does not contain the expected fields")
    colnames(ans)[1L] <- "assembly_accession"
    .check_assembly_summary(ans, genbank_or_refseq)
    assign(objname, ans, envir=.assembly_summary_cache)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### build_and_save_assembly_accessions_table()
###
### Use this utility to update assembly_accessions dataset located in
### GenomeInfoDb package (in /inst/extdata/).
###

.make_assembly_accessions_table_from_assembly_summaries <-
    function(genbank_assembly_summary, refseq_assembly_summary)
{
    genbank_accession0 <- genbank_assembly_summary[ , "assembly_accession"]
    genbank_rs_paired_asm0 <- genbank_assembly_summary[ , "gbrs_paired_asm"]
    refseq_accession0 <- refseq_assembly_summary[ , "assembly_accession"]
    refseq_gb_paired_asm0 <- refseq_assembly_summary[ , "gbrs_paired_asm"]

    df1 <- data.frame(genbank_accession=genbank_accession0,
                      refseq_accession=genbank_rs_paired_asm0,
                      stringsAsFactors=FALSE)
    df2 <- data.frame(genbank_accession=refseq_gb_paired_asm0,
                      refseq_accession=refseq_accession0,
                      stringsAsFactors=FALSE)
    ans <- merge.data.frame(df1, df2, all=TRUE)
    genbank_accession <- ans[ , "genbank_accession"]
    refseq_accession <- ans[ , "refseq_accession"]
    if (anyDuplicated(genbank_accession, incomparables=NA)
     || anyDuplicated(refseq_accession, incomparables=NA))
        stop("GenomeInfoDb internal error")

    m1 <- match(genbank_accession, genbank_accession0)
    m2 <- match(refseq_accession, refseq_accession0)
    paired_asm_comp1 <- genbank_assembly_summary[m1, "paired_asm_comp"]
    paired_asm_comp2 <- refseq_assembly_summary[m2, "paired_asm_comp"]
    is_na1 <- is.na(m1)
    is_na2 <- is.na(m2)
    stopifnot(all(is_na1 | is_na2 | paired_asm_comp1 == paired_asm_comp2))
    paired_asm_comp1[is_na1] <- paired_asm_comp2[is_na1]

    GCA_ftp_path <- genbank_assembly_summary[m1, "ftp_path"]
    GCF_ftp_path <- refseq_assembly_summary[m2, "ftp_path"]
    cbind(ans, paired_asm_comp=paired_asm_comp1,
               GCA_ftp_path=GCA_ftp_path,
               GCF_ftp_path=GCF_ftp_path,
               stringsAsFactors=FALSE)
}

### Saves "assembly_accessions.rda" in directory specified thru 'dir'.
build_and_save_assembly_accessions_table <- function(dir=".", quiet=FALSE)
{
    gb_assembly_summary <- fetch_assembly_summary("genbank", quiet=quiet)
    rs_assembly_summary <- fetch_assembly_summary("refseq", quiet=quiet)
    assembly_accessions <-
        .make_assembly_accessions_table_from_assembly_summaries(
            gb_assembly_summary, rs_assembly_summary)
    save(assembly_accessions, file=file.path(dir, "assembly_accessions.rda"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetch_assembly_report()
###

### 'assembly_accession' can be:
###   (a) a GenBank assembly accession (e.g. "GCA_000001405.15");
###   (b) an NCBI assembly name (e.g. "GRCh38").
.lookup_refseq_assembly_accession <- function(assembly_accession)
{
    genbank_assembly_summary <- fetch_assembly_summary("genbank")
    refseq_assembly_summary <- fetch_assembly_summary("refseq")

    .get_answers <- function(idx1, idx2) {
        ans2 <- refseq_assembly_summary$assembly_accession[idx2]
        if (all(is.na(idx1)))
            return(ans2)
        ans1 <- genbank_assembly_summary$gbrs_paired_asm[idx1]
        ans1 <- setdiff(ans1, "na")
        union(ans1, ans2)
    }

    if (.is_genbank_assembly_accession(assembly_accession)) {
        idx1 <- match(assembly_accession,
                      genbank_assembly_summary$assembly_accession)
        idx2 <- which(refseq_assembly_summary$gbrs_paired_asm ==
                      assembly_accession)
        ans <- .get_answers(idx1, idx2)
        if (length(ans) == 1L)
            return(ans)
        if (length(ans) >= 2L) {
            warning("more than one RefSeq assembly accession found for ",
                    "\"", assembly_accession, "\",\n  returning the 1st one")
            return(ans[1L])
        }
        return(NA_character_)
    }

    ## If 'assembly_accession' is not a RefSeq or GenBank assembly accession,
    ## then we assume it's an assembly name (e.g. "GRCh38", or "hg16", or
    ## "Pan_troglodytes-2.1.4").

    ## Exact match.
    idx1 <- which(genbank_assembly_summary$asm_name == assembly_accession)
    idx2 <- which(refseq_assembly_summary$asm_name == assembly_accession)
    ans <- .get_answers(idx1, idx2)
    if (length(ans) == 1L)
        return(ans)
    if (length(ans) >= 2L)
        stop("more than one RefSeq assembly accession found for ",
             "\"", assembly_accession, "\"")

    ## Fuzzy match.
    warning("No RefSeq assembly accession found for ",
            "\"", assembly_accession, "\".\n",
            "  Searching again using ",
            "\"", assembly_accession, "\" as a regular expression.")
    idx1 <- grep(assembly_accession, genbank_assembly_summary$asm_name,
                 ignore.case=TRUE)
    idx2 <- grep(assembly_accession, refseq_assembly_summary$asm_name,
                 ignore.case=TRUE)
    ans <- .get_answers(idx1, idx2)
    if (length(ans) == 1L)
        return(ans)
    if (length(ans) >= 2L) 
        stop("more than one RefSeq assembly accession found for regular ",
             "expression \"", assembly_accession, "\"")

    NA_character_
}

### 'assembly_accession' can be:
###   (a) a GenBank assembly accession (e.g. "GCA_000001405.15");
###   (b) a RefSeq assembly accession (e.g. "GCF_000001405.26");
###   (c) an NCBI assembly name (e.g. "GRCh38").
.normarg_assembly_accession <- function(assembly_accession)
{
    if (.is_genbank_assembly_accession(assembly_accession) ||
        .is_refseq_assembly_accession(assembly_accession))
        return(assembly_accession)
    accession0 <- assembly_accession
    assembly_accession <- .lookup_refseq_assembly_accession(accession0)
    if (is.na(assembly_accession))
        stop("cannot find a RefSeq assembly accession for \"",
             accession0, "\"")
    assembly_accession
}

.make_assembly_report_URL <- function(assembly_accession)
{
    assembly_accession <- .normarg_assembly_accession(assembly_accession)
    assembly_report_filename <- paste0(assembly_accession, ".assembly.txt")
    paste0(.NCBI_ASSEMBLY_REPORTS_URL, "All/", assembly_report_filename)
}

.fetch_assembly_report_from_URL <- function(url)
{
    ## R doesn't like seeing dashes ("-") or slashes ("/") in the colnames
    ## of a data frame. So we remove them from the official field names used
    ## by NCBI in the assembly report.
    colnames <- c("SequenceName", "SequenceRole", "AssignedMolecule",
                  "AssignedMoleculeLocationOrType", "GenBankAccn",
                  "Relationship", "RefSeqAccn", "AssemblyUnit" , 
                  "SequenceLength", "UCSCStyleName")
    read.table(url, sep="\t", col.names=colnames, stringsAsFactors=FALSE)
}

### See .normarg_assembly_accession() for how 'assembly_accession' can be
### specified. In addition, here 'assembly_accession' can be the URL to an
### assembly report (a.k.a. full sequence report). Examples of such URLs:
###   ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.26.assembly.txt
###   ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/GCA_000001405.15_GRCh38_assembly_report.txt
### Note that the 2 URls above both point to the assembly report for GRCh38,
### but the report at the 1st URL contains a little bit more information than
### the report at the 2nd URL.
fetch_assembly_report <- function(assembly_accession, AssemblyUnits=NULL)
{
    if (!isSingleString(assembly_accession))
        stop("'assembly_accession' must be a single string")
    if (grepl("://", assembly_accession, fixed=TRUE)) {
        report_url <- assembly_accession
    } else {
        report_url <- .make_assembly_report_URL(assembly_accession)
    }
    ans <- .fetch_assembly_report_from_URL(report_url)
    if (!is.null(AssemblyUnits)) {
        stopifnot(all(AssemblyUnits %in% ans$AssemblyUnit))
        idx <- which(ans$AssemblyUnit %in% AssemblyUnits)
        ans <- ans[idx, , drop=FALSE]
    }
    ans
}

