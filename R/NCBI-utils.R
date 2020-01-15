### =========================================================================
### Some low-level utilities to fetch data from NCBI
### -------------------------------------------------------------------------
###


.NCBI_ASSEMBLY_REPORTS_URL <-
    "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/"

.NCBI_ALL_ASSEMBLY_URL <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/all"
.GENBANK_ASSEMBLY_ACCESSION_PREFIX <- "GCA_"
.REFSEQ_ASSEMBLY_ACCESSION_PREFIX <- "GCF_"

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

.assembly_summary_cache <- new.env(parent=emptyenv())


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetch_assembly_summary()
###

.PAIRED_ASM_COMP_LEVELS <- c("identical", "different", NA)

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
    stopifnot(all(paired_asm_comp %in% .PAIRED_ASM_COMP_LEVELS))

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

### NOT exported.
fetch_assembly_summary <- function(genbank_or_refseq, quiet=FALSE)
{
    objname <- paste0("assembly_summary_", genbank_or_refseq)
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
### It will issue a warning that a small number of assemblies were dropped
### (3 on Feb 2017). It's OK to ignore if the number is small.
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
    is_na1 <- is.na(m1)
    is_na2 <- is.na(m2)

    ## Add "asm_name", "organism_name", and "taxid" columns.
    asm_name1 <- genbank_assembly_summary[m1, "asm_name"]
    organism_name1 <- genbank_assembly_summary[m1, "organism_name"]
    taxid1 <- genbank_assembly_summary[m1, "taxid"]
    ans1 <- data.frame(asm_name=asm_name1,
                       organism_name=organism_name1,
                       taxid=taxid1,
                       stringsAsFactors=FALSE)
    asm_name2 <- refseq_assembly_summary[m2, "asm_name"]
    organism_name2 <- refseq_assembly_summary[m2, "organism_name"]
    taxid2 <- refseq_assembly_summary[m2, "taxid"]
    ans2 <- data.frame(asm_name=asm_name2,
                       organism_name=organism_name2,
                       taxid=taxid2,
                       stringsAsFactors=FALSE)
    disagrement_idx <- which(!(is_na1 | is_na2 |
                               (asm_name1 == asm_name2 &
                                organism_name1 == organism_name2 &
                                taxid1 == taxid2)))
    if (length(disagrement_idx) != 0L) {
        warning(wmsg(length(disagrement_idx), " assemblies were dropped ",
                     "because the asm_name, organism_name, and/or taxid ",
                     "reported by GenBank and RefSeq do not agree"))
        ans <- ans[-disagrement_idx, , drop=FALSE]
        genbank_accession <- ans[ , "genbank_accession"]
        refseq_accession <- ans[ , "refseq_accession"]
        m1 <- m1[-disagrement_idx]
        m2 <- m2[-disagrement_idx]
        is_na1 <- is.na(m1)
        is_na2 <- is.na(m2)
        ans1 <- ans1[-disagrement_idx, , drop=FALSE]
        ans2 <- ans2[-disagrement_idx, , drop=FALSE]
    }
    ans1[is_na1, c("asm_name", "organism_name", "taxid")] <-
        ans2[is_na1, c("asm_name", "organism_name", "taxid")]
    ## Turning the "organism_name" column into a factor makes it bigger in
    ## memory and on disk. Not worth it.
    #ans1[ , "organism_name"] <- factor(ans1[ , "organism_name"],
    #                                   levels=unique(ans1[ , "organism_name"]))
    ans <- cbind(ans, ans1)

    ## Add "paired_asm_comp", "GCA_ftp_path", and "GCF_ftp_path" columns.
    paired_asm_comp1 <- genbank_assembly_summary[m1, "paired_asm_comp"]
    paired_asm_comp2 <- refseq_assembly_summary[m2, "paired_asm_comp"]
    stopifnot(all(is_na1 | is_na2 | paired_asm_comp1 == paired_asm_comp2))
    paired_asm_comp1[is_na1] <- paired_asm_comp2[is_na1]
    paired_asm_comp <- factor(paired_asm_comp1,
                              levels=.PAIRED_ASM_COMP_LEVELS)

    GCA_ftp_path <- genbank_assembly_summary[m1, "ftp_path"]
    GCF_ftp_path <- refseq_assembly_summary[m2, "ftp_path"]

    cbind(ans, paired_asm_comp=paired_asm_comp,
               GCA_ftp_path=GCA_ftp_path,
               GCF_ftp_path=GCF_ftp_path,
               stringsAsFactors=FALSE)
}

### NOT exported.
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
    objname <- "assembly_accessions"
    assembly_accessions <- try(get(objname, envir=.assembly_summary_cache,
                                   inherits=FALSE), silent=TRUE)
    if (is(assembly_accessions, "try-error")) {
        filename <- paste0(objname, ".rda")
        filepath <- system.file("extdata", filename, package="GenomeInfoDb")
        load(filepath, envir=.assembly_summary_cache)
        assembly_accessions <- get(objname, envir=.assembly_summary_cache,
                                   inherits=FALSE)
    }

    if (.is_genbank_assembly_accession(assembly_accession)) {
        idx <- match(assembly_accession,
                     assembly_accessions[ , "genbank_accession"])
        return(assembly_accessions[idx, "refseq_accession"])
    }

    ## If 'assembly_accession' is not a RefSeq or GenBank assembly accession,
    ## then we assume it's an assembly name (e.g. "GRCh38", or "hg16", or
    ## "Pan_troglodytes-2.1.4").

    ## Exact match.
    idx <- which(assembly_accessions[ , "asm_name"] == assembly_accession)
    if (length(idx) == 1L)
        return(assembly_accessions[idx , "refseq_accession"])
    if (length(idx) >= 2L)
        stop("more than one RefSeq assembly accession found for ",
             "\"", assembly_accession, "\"")

    ## Fuzzy match.
    warning("No RefSeq assembly accession found for ",
            "\"", assembly_accession, "\".\n",
            "  Searching again using ",
            "\"", assembly_accession, "\" as a regular expression.")
    idx <- grep(assembly_accession, assembly_accessions[ , "asm_name"],
                ignore.case=TRUE)
    if (length(idx) == 1L)
        return(assembly_accessions[idx, "refseq_accession"])
    if (length(idx) >= 2L)
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

### Returns https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_report.txt for GCA_000001405.15
.make_assembly_report_URL <- function(assembly_accession)
{
    assembly_accession <- .normarg_assembly_accession(assembly_accession)
    prefix <- substr(assembly_accession,
                     1L, nchar(.GENBANK_ASSEMBLY_ACCESSION_PREFIX))
    if (prefix == .GENBANK_ASSEMBLY_ACCESSION_PREFIX) {
        url <- "GCA"
    } else if (prefix == .REFSEQ_ASSEMBLY_ACCESSION_PREFIX) {
        url <- "GCF"
    } else {
        stop(wmsg("don't know where to find assembly report for ",
                  assembly_accession))
    }
    parts_end <- nchar(prefix) + (1:3) * 3L
    parts <- substring(assembly_accession, parts_end - 2L, parts_end)
    url <- paste0(.NCBI_ALL_ASSEMBLY_URL, "/", url, "/",
                  paste0(parts, collapse="/"), "/")
    listing <- list_ftp_dir(url)
    idx <- which(paste0(assembly_accession, "_") ==
                 substr(listing, 1L, nchar(assembly_accession)+1L))
    if (length(idx) == 0L)
        stop(wmsg("don't know where to find assembly report for ",
                  assembly_accession))
    if (length(idx) > 1L)
        stop(wmsg("more than one assembly report found for ",
                  assembly_accession))
    part4 <- listing[[idx]]
    paste0(url, part4, "/", part4, "_assembly_report.txt")
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
    read.table(url, sep="\t", col.names=colnames, na.strings=c("NA", "na"),
               stringsAsFactors=FALSE)
}

### NOT exported.
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load and access db of registered NCBI genomes
###

.NCBI_genome2accession <- new.env(parent=emptyenv())
.NCBI_accession2assembly <- new.env(parent=emptyenv())

.load_registered_NCBI_genome <- function(file_path)
{
    filename <- basename(file_path)
    if (substr(filename, nchar(filename)-1L, nchar(filename)) != ".R")
        stop(wmsg("name of genome registration file '", filename, "' must ",
                  "have extension .R"))
    if (grepl(" ", filename, fixed=TRUE))
        stop(wmsg("name of genome registration file '", filename, "' must ",
                  "not contain spaces"))

    ## Placeholders. Will actually get defined when we source the
    ## assembly files.
    ORGANISM <- NULL    # Expected to be a single string.
    ASSEMBLIES <- NULL  # Expected to be a list with one list element per
                        # assembly.
    source(file_path, local=TRUE)

    stop_if <- function(notok, ...) {
        if (notok)
            stop("Error in NCBI genome registration file '", filename,
                 "':\n  ", wmsg(...))
    }

    ## Check ORGANISM.
    stop_if(is.null(ORGANISM), "'ORGANISM' must be defined")
    stop_if(!isSingleString(ORGANISM), "'ORGANISM' must be a single string")
    stop_if(grepl("_", ORGANISM, fixed=TRUE),
            "underscores are not allowed in 'ORGANISM' (use spaces instead)")
    target <- chartr("_", " ", substr(filename, 1L, nchar(filename)-2L))
    stop_if(!identical(target, ORGANISM),
            "'ORGANISM' must match filename ",
            "(with spaces in place of underscores)")

    ## Check ASSEMBLIES.
    stop_if(is.null(ASSEMBLIES), "'ASSEMBLIES' must be defined")
    stop_if(!is.list(ASSEMBLIES), "'ASSEMBLIES' must be a list of lists")

    required_fields <- c("genome", "assembly_accession", "date", "circ_seqs")
    for (i in seq_along(ASSEMBLIES)) {
        assembly <- ASSEMBLIES[[i]]

        stop_if(!is.list(assembly),
                "'ASSEMBLIES[[", i, "]]' must be a named list")
        assembly_fields <- names(assembly)
        stop_if(!is.character(assembly_fields),
                "'ASSEMBLIES[[", i, "]]' must be a named list")
        stop_if(anyNA(assembly_fields) ||
                !all(nzchar(assembly_fields)) ||
                anyDuplicated(assembly_fields),
                "the names on 'ASSEMBLIES[[", i, "]]' must ",
                "not contain NAs, empty strings, or duplicates")
        stop_if(!all(required_fields %in% names(assembly)),
                "'ASSEMBLIES[[", i, "]]' must have fields: ",
                paste(paste0("\"", required_fields, "\""), collapse=", "))

        ## Check "genome" field (required).
        genome <- assembly$genome
        stop_if(!isSingleString(genome) || genome == "",
                "\"genome\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "be a single non-empty string")

        ## Check "accession" field (required).
        accession <- assembly$assembly_accession
        stop_if(!isSingleString(accession) || accession == "",
                "\"accession\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "be a single non-empty string")
        stop_if(!is.null(.NCBI_accession2assembly[[accession]]),
                "assembly accession ", accession, " used more than once")

        ## Check "date" field (required).
        date <- assembly$date
        stop_if(!isSingleString(date) || date == "",
                "\"date\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "be a single non-empty string")

        ## Check "circ_seqs" field (required).
        circ_seqs <- assembly$circ_seqs
        stop_if(!is.character(circ_seqs),
                "\"circ_seqs\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "be a character vector")
        stop_if(anyNA(circ_seqs) ||
                !all(nzchar(circ_seqs)) ||
                anyDuplicated(circ_seqs),
                "\"circ_seqs\" field in 'ASSEMBLIES[[", i, "]]' must ",
                "not contain NAs, empty strings, or duplicates")

        ## Check optional fields.
        infraspecific_name <- assembly$infraspecific_name
        if (!is.null(infraspecific_name)) {
            stop_if(!isSingleString(infraspecific_name) ||
                    infraspecific_name == "" ||
                    !isSingleString(names(infraspecific_name)),
              "\"infraspecific_name\" field in 'ASSEMBLIES[[", i, "]]' must ",
              "be a single non-empty named string")
        }

        genome <- tolower(genome)
        .NCBI_genome2accession[[genome]] <-
            c(.NCBI_genome2accession[[genome]], accession)
        assembly$organism <- ORGANISM
        .NCBI_accession2assembly[[accession]] <- assembly
    }
}

.load_registered_NCBI_genomes <- function()
{
    dir_path <- system.file("registered_genomes", "NCBI",
                             package="GenomeInfoDb")
    file_paths <- list.files(dir_path, pattern="\\.R$", full.names=TRUE)
    for (file_path in file_paths)
        .load_registered_NCBI_genome(file_path)
}

registered_NCBI_genomes <- function()
{
    if (length(.NCBI_accession2assembly) == 0L)
        .load_registered_NCBI_genomes()
    assemblies <- unname(as.list(.NCBI_accession2assembly, all.names=TRUE))

    colnames <- c("organism", "genome", "assembly_accession",
                  "infraspecific_name", "date", "circ_seqs")
    make_col <- function(colname) {
        col0 <- lapply(assemblies, `[[`, colname)
        if (colname == "circ_seqs")
            return(CharacterList(col0))
        if (colname == "infraspecific_name") {
            stopifnot(all(lengths(col0) <= 1L))
            col <- rep.int(NA_character_, length(col0))
            idx1 <- which(lengths(col0) == 1L)
            if (length(idx1) != 0L) {
                col1 <- unlist(col0)
                stopifnot(is.character(col1))
                col1_names <- names(col1)
                stopifnot(is.character(col1_names),
                          !anyNA(col1_names),
                          all(nzchar(col1_names)))
                col[idx1] <- paste(col1_names, col1, sep=":")
            }
            return(factor(col))
        }
        stopifnot(all(lengths(col0) == 1L))
        col <- as.character(unlist(col0, use.names=FALSE))
        if (colname == "organism")
            col <- factor(col)  # order of levels will dictate order
                                # of rows in final DataFrame
        col
    }

    listData <- lapply(setNames(colnames, colnames), make_col)
    ans <- S4Vectors:::new_DataFrame(listData, nrows=length(assemblies))
    oo <- order(ans$organism, ans$infraspecific_name, ans$date)
    as.data.frame(ans[oo, , drop=FALSE])
}

lookup_NCBI_genome2accession <- function(genome)
{
    stopifnot(isSingleString(genome))
    if (length(.NCBI_genome2accession) == 0L)
        .load_registered_NCBI_genomes()
    .NCBI_genome2accession[[tolower(genome)]]
}

lookup_NCBI_accession2assembly <- function(accession)
{
    stopifnot(isSingleString(accession))
    if (length(.NCBI_accession2assembly) == 0L)
        .load_registered_NCBI_genomes()
    .NCBI_accession2assembly[[accession]]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_chrom_info_from_NCBI()
###

.format_assembly_report <- function(assembly_report, circ_seqs=NULL)
{
    drop_columns <- c("AssignedMolecule", "AssignedMoleculeLocationOrType")
    ans <- assembly_report[ , !(colnames(assembly_report) %in% drop_columns)]
    ans[ , "SequenceName"] <- as.character(ans[ , "SequenceName"])
    SequenceRole_levels <- c("assembled-molecule",
                             "alt-scaffold",
                             "unlocalized-scaffold",
                             "unplaced-scaffold",
                             "pseudo-scaffold")
    sequence_role <- factor(ans[ , "SequenceRole"], levels=SequenceRole_levels)
    ans[ , "SequenceRole"] <- sequence_role
    oo <- order(as.integer(sequence_role))
    ans <- S4Vectors:::extract_data_frame_rows(ans, oo)
    Relationship_levels <- c("=", "<>")
    ans[ , "Relationship"] <- factor(ans[ , "Relationship"],
                                     levels=Relationship_levels)
    ans[ , "AssemblyUnit"] <- factor(ans[ , "AssemblyUnit"])
    na_idx <- which(ans[ , "UCSCStyleName"] %in% "na")
    ans[na_idx , "UCSCStyleName"] <- NA_character_
    circular <- make_circ_flags_from_circ_seqs(ans[ , "SequenceName"],
                                               circ_seqs=circ_seqs)
    stopifnot(all(ans[which(circular), "SequenceRole"] %in%
                  "assembled-molecule"))
    ans$circular <- circular
    ans
}

.NCBI_cached_chrom_info <- new.env(parent=emptyenv())

.get_NCBI_chrom_info_from_accession <- function(accession, circ_seqs=NULL,
    assembled.molecules.only=FALSE,
    recache=FALSE)
{
    ans <- .NCBI_cached_chrom_info[[accession]]
    if (is.null(ans) || recache) {
        assembly_report <- fetch_assembly_report(accession)
        ans <- .format_assembly_report(assembly_report, circ_seqs=circ_seqs)
        .NCBI_cached_chrom_info[[accession]] <- ans
    }
    if (assembled.molecules.only) {
        keep_idx <- which(ans[ , "SequenceRole"] %in% "assembled-molecule")
        ans <- S4Vectors:::extract_data_frame_rows(ans, keep_idx)
    }
    ans
}

### Return an 8-column data.frame with columns "SequenceName" (character),
### "SequenceRole" (factor),  "GenBankAccn" (character), "Relationship"
### (factor), "RefSeqAccn" (character), "AssemblyUnit" (factor),
### "SequenceLength" (integer), "UCSCStyleName" (character).
get_chrom_info_from_NCBI <- function(genome,
    assembled.molecules.only=FALSE,
    recache=FALSE)
{
    if (!isSingleString(genome))
        stop(wmsg("'genome' must be a single string"))
    if (!isTRUEorFALSE(assembled.molecules.only))
        stop(wmsg("'assembled.molecules.only' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))

    ## First see if the user supplied the assembly accession of a registered
    ## genome assembly instead of its name.
    assembly <- lookup_NCBI_accession2assembly(genome)
    if (!is.null(assembly)) {
        ## Yes s/he did.
        accession <- genome
        circ_seqs <- assembly$circ_seqs
    } else {
        ## No s/he didn't.
        ## Now see if s/he supplied the name of a registered genome assembly.
        accession <- lookup_NCBI_genome2accession(genome)
        if (!is.null(accession)) {
            ## Yes s/he did.
            if (length(accession) > 1L) {
                in1string <- paste0(accession, collapse=", ")
                warning(wmsg("Genome ", genome, " is mapped to more ",
                             "then one assembly (", in1string, "). ",
                             "The first one was selected."))
                accession <- accession[[1L]]
            }
            assembly <- lookup_NCBI_accession2assembly(accession)
            circ_seqs <- assembly$circ_seqs
        } else {
            ## No s/he didn't.
            ## So now we just assume that 'genome' is an assembly accession
            ## (an unregistered one).
            accession <- genome
            circ_seqs <- NULL
        }
    }
    .get_NCBI_chrom_info_from_accession(accession,
            circ_seqs=circ_seqs,
            assembled.molecules.only=assembled.molecules.only,
            recache=recache)
}

