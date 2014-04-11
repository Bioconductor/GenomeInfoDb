### =========================================================================
### Some low-level (non exported) utility functions.
### -------------------------------------------------------------------------


.GENBANK_ASSEMBLY_ID_PREFIX <- "GCA_"
.REFSEQ_ASSEMBLY_ID_PREFIX <- "GCF_"

.NCBI_ASSEMBLY_REPORTS_URL <-
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/"

.GENBANK_GENOMES_URL <-
    "ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/"

.isSingleString <- function(x)
{
    is.character(x) && length(x) == 1L && !is.na(x)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetch_assembly_report()
###

.is_genbank_assembly_id <- function(x)
{
    ## We use %in% instead of == to be NA-proof.
    substr(x, 1L, nchar(.GENBANK_ASSEMBLY_ID_PREFIX)) %in%
        .GENBANK_ASSEMBLY_ID_PREFIX
}

.is_refseq_assembly_id <- function(x)
{
    ## We use %in% instead of == to be NA-proof.
    substr(x, 1L, nchar(.REFSEQ_ASSEMBLY_ID_PREFIX)) %in%
        .REFSEQ_ASSEMBLY_ID_PREFIX
}

### Performs some quick sanity checks on the assembly summary.
.check_assembly_summary <- function(assembly_summary, genbank_or_refseq)
{
    ## Check "assembly_id" prefixes.
    if (genbank_or_refseq == "genbank") {
        is_valid_id <- .is_genbank_assembly_id
    } else {
        is_valid_id <- .is_refseq_assembly_id
    }
    stopifnot(all(is_valid_id(assembly_summary$assembly_id)))

    ## Check "assembly_id" uniqueness.
    stopifnot(anyDuplicated(assembly_summary$assembly_id) == 0L)

    ## Check "paired_asm_comp" levels.
    paired_asm_comp_levels <- c("identical", "different", "na")
    stopifnot(all(assembly_summary$paired_asm_comp %in%
                  paired_asm_comp_levels))

    ## Check gbrs_paired_asm/paired_asm_comp consistency.
    gbrs_paired_asm <- assembly_summary$gbrs_paired_asm
    idx1 <- which(gbrs_paired_asm %in% "na")
    idx2 <- which(assembly_summary$paired_asm_comp %in% "na")
    stopifnot(identical(idx1, idx2))
    if (length(idx1) != 0L)
        gbrs_paired_asm <- gbrs_paired_asm[-idx1]

    ## Check "gbrs_paired_asm" prefixes.
    if (genbank_or_refseq == "genbank") {
        is_valid_id <- .is_refseq_assembly_id
    } else {
        is_valid_id <- .is_genbank_assembly_id
    }
    stopifnot(all(is_valid_id(gbrs_paired_asm)))
}

.assembly_summary_cache <- new.env(parent=emptyenv())

.fetch_assembly_summary <- function(genbank_or_refseq)
{
    objname <- paste0("assembly_summary_", genbank_or_refseq)
    ans <- try(get(objname, envir=.assembly_summary_cache, inherits=FALSE),
               silent=TRUE)
    if (!is(ans, "try-error"))
        return(ans)
    assembly_summary_filename <- paste0(objname, ".txt")
    url <- paste0(.NCBI_ASSEMBLY_REPORTS_URL, assembly_summary_filename)
    destfile <- tempfile()
    download.file(url, destfile, quiet=TRUE)
    ans <- read.table(destfile, header=TRUE, sep="\t", quote="",
                                comment.char="", stringsAsFactors=FALSE)
    target_colnames <- c("X..assembly_id", "bioproject", "biosample",
                         "wgs_master", "representative_status", "taxid",
                         "species_taxid", "organism_name",
                         "infraspecific_name", "isolate", "version_status",
                         "assembly_level", "release_type", "genome_rep",
                         "seq_rel_date", "asm_name", "submitter",
                         "gbrs_paired_asm", "paired_asm_comp")
    if (!identical(target_colnames, colnames(ans)))
        stop(url, " does not contain the expected fields")
    colnames(ans)[1L] <- "assembly_id"
    .check_assembly_summary(ans, genbank_or_refseq)
    assign(objname, ans, envir=.assembly_summary_cache)
    ans
}

### 'assembly' can be:
###   (a) a RefSeq Assembly ID (e.g. "GCF_000001405.26"), in which case it's
###       returned as is;
###   (b) a GenBank Assembly ID (e.g. "GCA_000001405.15");
###   (c) an NCBI assembly name (e.g. "GRCh38").
lookup_refseq_assembly_id <- function(assembly)
{
    if (!.isSingleString(assembly))
        stop("'assembly' must be a single string")
    if (.is_refseq_assembly_id(assembly))
        return(assembly)

    assembly_summary_genbank <- .fetch_assembly_summary("genbank")
    assembly_summary_refseq <- .fetch_assembly_summary("refseq")

    .get_answers <- function(idx1, idx2) {
        ans1 <- assembly_summary_genbank$gbrs_paired_asm[idx1]
        ans1 <- setdiff(ans1, "na")
        ans2 <- assembly_summary_refseq$assembly_id[idx2]
        union(ans1, ans2)
    }

    if (.is_genbank_assembly_id(assembly)) {
        idx1 <- match(assembly, assembly_summary_genbank$assembly_id)
        idx2 <- which(assembly_summary_refseq$gbrs_paired_asm == assembly)
        ans <- .get_answers(idx1, idx2)
        if (length(ans) == 1L)
            return(ans)
        if (length(ans) >= 2L) {
            warning("more than one RefSeq assembly id found for ",
                    "\"", assembly, "\",\n  returning the 1st one")
            return(ans[1L])
        }
        return(NA_character_)
    }

    ## If 'assembly' is not a RefSeq Assembly ID or a GenBank Assembly ID,
    ## then we assume it's an assembly name (e.g. "GRCh38" or
    ## "Pan_troglodytes-2.1.4").

    ## Exact match.
    idx1 <- which(assembly_summary_genbank$asm_name == assembly)
    idx2 <- which(assembly_summary_refseq$asm_name == assembly)
    ans <- .get_answers(idx1, idx2)
    if (length(ans) == 1L)
        return(ans)
    if (length(ans) >= 2L)
        stop("more than one RefSeq assembly id found for ",
             "\"", assembly, "\"")

    ## Fuzzy match.
    warning("No RefSeq assembly id found for ",
            "\"", assembly, "\".\n",
            "  Searching again using ",
            "\"", assembly, "\" as a regular expression.")
    idx1 <- grep(assembly, assembly_summary_genbank$asm_name, ignore.case=TRUE)
    idx2 <- grep(assembly, assembly_summary_refseq$asm_name, ignore.case=TRUE)
    ans <- .get_answers(idx1, idx2)
    if (length(ans) == 1L)
        return(ans)
    if (length(ans) >= 2L) 
        stop("more than one RefSeq assembly id found for regular ",
             "expression \"", assembly, "\"")

    NA_character_
}

### See lookup_refseq_assembly_id() for how 'assembly' can be specified.
.make_assembly_report_URL <- function(assembly)
{
    assembly_id <- lookup_refseq_assembly_id(assembly)
    assembly_report_filename <- paste0(assembly_id, ".assembly.txt")
    paste0(.NCBI_ASSEMBLY_REPORTS_URL, "All/", assembly_report_filename)
}

.fetch_assembly_report_from_URL <- function(url)
{
    ## R doesn't like seeing dashes ("-") or slashes ("/") in the colnames
    ## of a data frame. So we remove them from the official field names used
    ## by NCBI in the assembly report.
    colnames <- c("SequenceName", "SequenceRole", "AssignedMolecule",
                  "AssignedMoleculeLocationOrType", "GenBankAccn",
                  "Relationship", "RefSeqAccn", "AssemblyUnit")
    read.table(url, sep="\t", col.names=colnames, stringsAsFactors=FALSE)
}

### See lookup_refseq_assembly_id() for how 'assembly' can be specified.
### In addition, here 'assembly' can be the URL to an assembly report (a.k.a.
### full sequence report). Examples of such URLs:
###   ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.26.assembly.txt
###   ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/GCA_000001405.15_GRCh38_assembly_report.txt
### Note that the 2 URls above both point to the assembly report for GRCh38,
### but the report at the 1st URL contains a little bit more information than
### the report at the 2nd URL.
fetch_assembly_report <- function(assembly, AssemblyUnits=NULL)
{
    if (!.isSingleString(assembly))
        stop("'assembly' must be a single string")
    if (!grepl("://", assembly))
        assembly <- .make_assembly_report_URL(assembly)
    ans <- .fetch_assembly_report_from_URL(assembly)
    if (!is.null(AssemblyUnits)) {
        stopifnot(all(AssemblyUnits %in% ans$AssemblyUnit))
        idx <- which(ans$AssemblyUnit %in% AssemblyUnits)
        ans <- ans[idx, , drop=FALSE]
    }
    ans
}

