### =========================================================================
### Some low-level utilities to fetch data from Ensembl
### -------------------------------------------------------------------------
###
### Unless stated otherwise, nothing in this file is exported.
###
### List of Ensembl public MySQL servers / ports
###   https://www.ensembl.org/info/data/mysql.html
### Ensembl Core Schema:
###   https://www.ensembl.org/info/docs/api/core/core_schema.html
### Full schema:
###   ftp://ftp.ensembl.org/pub/ensembl/sql/table.sql
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ensembl db schema (only for the tables we care about)
###

### Fundamental Tables

### 'attrib_type' table.
### nrow = 331 (in dataset hsapiens_gene_ensembl, Ensembl 99)
.ENSEMBLDB_ATTRIB_TYPE_COLUMNS <- c(
    "attrib_type_id",         # primary key
    "code",
    "name",
    "description"
)

### Assembly Tables

### 'seq_region' table.
### nrow = 268443 (in dataset hsapiens_gene_ensembl, Ensembl 99)
.ENSEMBLDB_SEQ_REGION_COLUMNS <- c(
    "seq_region_id",          # primary key
    "name",
    "coord_system_id",        # ==> coord_system.coord_system_id
    "length"
)

### 'coord_system' table.
### nrow = 9 (in dataset hsapiens_gene_ensembl, Ensembl 99)
.ENSEMBLDB_COORD_SYSTEM_COLUMNS <- c(
    "coord_system_id",        # primary key
    "species_id",
    "name",
    "version",
    "rank",
    "attrib"
)

### 'seq_region_attrib' table.
### nrow = 5927 (in dataset hsapiens_gene_ensembl, Ensembl 99)
.ENSEMBLDB_SEQ_REGION_ATTRIB_COLUMNS <- c(
    "seq_region_id",          # ==> seq_region.seq_region_id
    "attrib_type_id",         # ==> attrib_type.attrib_type_id
    "value"
)

### 'seq_region_synonym' table.
### nrow = 2360 (in dataset hsapiens_gene_ensembl, Ensembl 99)
.ENSEMBLDB_SEQ_REGION_SYNONYM_COLUMNS <- c(
    "seq_region_synonym_id",  # primary key
    "seq_region_id",          # ==> seq_region.seq_region_id
    "synonym",
    "external_db_id"          # ==> external_db.external_db_id
)

### External References

### 'external_db' table.
### nrow = 446 (in dataset hsapiens_gene_ensembl, Ensembl 99)
.ENSEMBLDB_EXTERNAL_DB_COLUMNS <- c(
    "external_db_id",         # primary key
    "db_name",
    "db_release",
    "status",
    "priority",
    "db_display_name",
    "type",
    "secondary_db_name",
    "secondary_db_table",
    "description"
)

.ENSEMBLDB_COLUMNS <- list(
    attrib_type=.ENSEMBLDB_ATTRIB_TYPE_COLUMNS,
    seq_region=.ENSEMBLDB_SEQ_REGION_COLUMNS,
    coord_system=.ENSEMBLDB_COORD_SYSTEM_COLUMNS,
    seq_region_attrib=.ENSEMBLDB_SEQ_REGION_ATTRIB_COLUMNS,
    seq_region_synonym=.ENSEMBLDB_SEQ_REGION_SYNONYM_COLUMNS,
    external_db=.ENSEMBLDB_EXTERNAL_DB_COLUMNS
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Fetch data from the Ensembl FTP server (via RCurl)
###
### The utililities in this section use RCurl or just utils::download.file()
### for getting stuff directly from the Ensembl FTP server. They can access
### stuff that is not available thru biomaRt like for example the lengths of
### the sequences in the reference genome associated with a particular dataset
### and Ensembl release (e.g. for dataset "hsapiens_gene_ensembl" and Ensembl
### release "64").
###
### Note that querying the Ensembl MySQL server (via RMariaDB) would probably
### be a better way to access this stuff.
### Update (Feb 10, 2020): Some preliminary testing indicates that using
### RMariaDB to fetch full tables is actually significantly slower.
### For example, to fetch table "seq_region" from db "homo_sapiens_core_99_38"
### (the table has 268443 rows and 4 columns):
###   o The utils::download.file() + utils::read.table() method takes about
###     7 sec to download and parse the MySQL dump located at:
###         ftp://ftp.ensembl.org/pub/current_mysql/homo_sapiens_core_99_38/
###   o The RMariaDB method takes about 22 sec to retrieve the table from the
###     MySQL server at ensembldb.ensembl.org.
###

.ENSEMBL.PUB_FTP_URL <- "ftp://ftp.ensembl.org/pub/"
.ENSEMBLGRCh37.PUB_FTP_URL <- "ftp://ftp.ensembl.org/pub/grch37/"
.ENSEMBLGENOMES.PUB_FTP_URL <- "ftp://ftp.ensemblgenomes.org/pub/"

### Uses the utils::download.file() + utils::read.table() method.
fetch_table_from_Ensembl_ftp <- function(core_url, table, full.colnames=FALSE,
                                         nrows=-1L)
{
    columns <- .ENSEMBLDB_COLUMNS[[table]]
    if (is.null(columns)) {
        warning(wmsg("unknown table: ", table, " (download might fail or ",
                     "the returned data frame will have automatic colnames)"))
    } else if (full.colnames) {
        columns <- paste(table, columns, sep=".")
    }
    url <- paste0(core_url, table, ".txt.gz")
    ans <- fetch_table_from_url(url, colnames=columns, nrows=nrows)
    if (full.colnames && is.null(columns))
        colnames(ans) <- paste(table, colnames(ans), sep=".")
    ans
}

### 'kingdom' must be NA or one of the EnsemblGenomes marts i.e. "bacteria",
### "fungi", "metazoa", "plants", or "protists".
get_url_to_Ensembl_ftp_mysql <- function(
    release=NA,
    use.grch37=FALSE, kingdom=NA)
{
    if (is.na(kingdom)) {
        if (is.na(release)) {
            if (use.grch37) {
                pub_subdir <- "current/mysql"
            } else {
                pub_subdir <- "current_mysql"
            }
        } else {
            pub_subdir <- paste0("release-", release, "/mysql")
        }
        if (use.grch37)
            pub_ftp_url <- .ENSEMBLGRCh37.PUB_FTP_URL
        else
            pub_ftp_url <- .ENSEMBL.PUB_FTP_URL
    } else {
        pub_ftp_url <- paste0(.ENSEMBLGENOMES.PUB_FTP_URL, kingdom, "/")
        if (is.na(release)) {
            pub_subdir <- "current"
        } else {
            pub_subdir <- paste0("release-", release)
        }
        pub_subdir <- paste0(pub_subdir, "/mysql")
    }
    paste0(pub_ftp_url, pub_subdir, "/")
}

get_url_to_Ensembl_ftp_gtf <- function(release=NA)
{
    if (is.na(release))
        pub_subdir <- "current_gtf"
    else
        pub_subdir <- paste0("release-", release, "/gtf")
    paste0(.ENSEMBL.PUB_FTP_URL, pub_subdir, "/")
}

### 'kingdom' must be NA or one of the EnsemblGenomes marts i.e. "bacteria",
### "fungi", "metazoa", "plants", or "protists".
list_Ensembl_ftp_mysql_core_dirs <- function(
    release=NA,
    use.grch37=FALSE, kingdom=NA, url=NA)
{
    if (is.na(url))
        url <- get_url_to_Ensembl_ftp_mysql(release, use.grch37, kingdom)
    core_dirs <- list_ftp_dir(url, subdirs.only=TRUE)
    pattern <- "_core_"
    if (!is.na(release))
        pattern <- paste0(pattern, release, "_")
    core_dirs[grep(pattern, core_dirs, fixed=TRUE)]
}

get_Ensembl_ftp_mysql_core_dir <- function(
    dataset, release=NA,
    use.grch37=FALSE, kingdom=NA, url=NA)
{
    if (is.na(url))
        url <- get_url_to_Ensembl_ftp_mysql(release, use.grch37, kingdom)
    core_dirs <- list_Ensembl_ftp_mysql_core_dirs(release=release,
                                                  use.grch37=use.grch37,
                                                  kingdom=kingdom,
                                                  url=url)
    trimmed_core_dirs <- sub("_core_.*$", "", core_dirs)
    shortnames <- sub("^(.)[^_]*_", "\\1", trimmed_core_dirs)
    if (dataset == "mfuro_gene_ensembl") {
        shortname0 <- "mputorius_furo"
    } else {
        shortname0 <- strsplit(dataset, "_", fixed=TRUE)[[1L]][1L]
    }
    core_dir <- core_dirs[shortnames == shortname0]
    if (length(core_dir) != 1L)
        stop("found 0 or more than 1 subdir for \"", dataset,
             "\" dataset at ", url)
    core_dir
}

### 'kingdom' must be NA or one of the EnsemblGenomes marts i.e. "bacteria",
### "fungi", "metazoa", "plants", or "protists".
### Return URL of Ensemble Core DB (FTP access).
get_url_to_Ensembl_ftp_mysql_core <- function(
    dataset, release=NA,
    use.grch37=FALSE, kingdom=NA, url=NA)
{
    if (is.na(url))
        url <- get_url_to_Ensembl_ftp_mysql(release, use.grch37, kingdom)
    core_dir <- get_Ensembl_ftp_mysql_core_dir(dataset, release=release,
                                               use.grch37=use.grch37,
                                               kingdom=kingdom,
                                               url=url)
    paste0(url, core_dir, "/")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetch_default_coord_systems_from_Ensembl_ftp()
###

fetch_default_coord_systems_from_Ensembl_ftp <- function(core_url)
{
    coord_system <- fetch_table_from_Ensembl_ftp(core_url, "coord_system")
    ## Drop rows that do not have the default_version attrib.
    keep_idx <- grep("default_version", coord_system[ , "attrib"], fixed=TRUE)
    S4Vectors:::extract_data_frame_rows(coord_system, keep_idx)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetch_seq_regions_from_Ensembl_ftp()
###
### This is the workhorse behind getChromInfoFromEnsembl().
###

### Retrieves all synonyms for the supplied sequence ids.
### This is done via tables "seq_region_synonym" and "external_db".
### Return a list parallel to 'seq_region_ids' where each list element is a
### named character vector containing all the synonyms for the corresponding
### sequence id. The names on the character vector indicate the origin of
### the synonym i.e. the external db where it's used (e.g. "INSDC",
### "RefSeq_genomic", "UCSC", "ensembl_internal_synonym", etc..)
.fetch_synonyms_from_Ensembl_ftp <- function(core_url, seq_region_ids)
{
    stopifnot(is_primary_key(seq_region_ids))

    all_synonyms <- fetch_table_from_Ensembl_ftp(core_url, "seq_region_synonym")
    keep_idx <- which(all_synonyms[ , "seq_region_id"] %in% seq_region_ids)
    synonyms <- S4Vectors:::extract_data_frame_rows(all_synonyms, keep_idx)

    external_dbs <- fetch_table_from_Ensembl_ftp(core_url, "external_db")
    synonyms <- join_dfs(synonyms, external_dbs,
                         "external_db_id", "external_db_id",
                         keep.Rcol=TRUE)  # we'll drop it below

    f <- factor(synonyms[ , "seq_region_id"], levels=seq_region_ids)
    unlisted_ans <- setNames(synonyms[ , "synonym"], synonyms[ , "db_name"])
    split(unlisted_ans, f)
}

.attrib_type_codes_to_ids_from_Ensembl_ftp <- function(core_url, codes)
{
    if (!is.character(codes))
        stop(wmsg("'codes' must be a character vector"))
    if (anyDuplicated(codes))
        stop(wmsg("'codes' cannot contain duplicates"))
    ## Get the first 99 rows of 'attrib_type' table.
    ## The reason we don't read in the entire table is because some rows at
    ## the bottom of the table (rows with attrib_type_id >= 416, these rows
    ## were added in Ensembl 75) contain embedded EOL characters that break
    ## read.table(). Since in practice we only need to lookup the
    ## attrib_type_id associated with the "toplevel" and "non_ref" codes,
    ## and since these codes are generally found at the beginning of the
    ## file (AFAIK "toplevel" has always been assigned the attrib_type_id
    ## of 6 and the lines in the file seem to always be ordered by
    ## attrib_type_id), reading in the first 99 rows should be way enough
    ## to get what we need.
    attrib_type <- fetch_table_from_Ensembl_ftp(core_url, "attrib_type",
                                                nrows=99L)
    m <- solid_match(codes, attrib_type[ , "code"],
                     x_what="supplied code",
                     table_what="\"attrib_type.code\" value")
    bad_idx <- which(is.na(m))
    if (length(bad_idx) != 0L) {
        in1string <- paste0(codes[bad_idx], collapse=", ")
        stop(wmsg("table \"attrib_type\" has no entry for: ", in1string))
    }
    setNames(attrib_type[m, "attrib_type_id"], codes)
}

.external_db_names_to_ids_from_Ensembl_ftp <- function(core_url, db_names)
{
    if (!is.character(db_names))
        stop(wmsg("'db_names' must be a character vector"))
    if (anyDuplicated(db_names))
        stop(wmsg("'db_names' cannot contain duplicates"))
    external_db <- fetch_table_from_Ensembl_ftp(core_url, "external_db")
    m <- solid_match(db_names, external_db[ , "db_name"],
                     x_what="supplied db_name",
                     table_what="\"external_db.db_name\" value")
    bad_idx <- which(is.na(m))
    if (length(bad_idx) != 0L) {
        in1string <- paste0(db_names[bad_idx], collapse=", ")
        stop(wmsg("table \"external_db\" has no entry for: ", in1string))
    }
    setNames(external_db[m, "external_db_id"], db_names)
}

### Retrieves attribs "toplevel" and "non_ref" for the supplied sequence ids.
### This is done via tables "seq_region_attrib" and "attrib_type".
.fetch_attribs_from_Ensembl_ftp <- function(core_url, seq_region_ids,
                                            toplevel=FALSE, non_ref=FALSE)
{
    all_attribs <- fetch_table_from_Ensembl_ftp(core_url, "seq_region_attrib")
    attrib_type_id <- all_attribs[ , "attrib_type_id"]
    codes <- c("toplevel", "non_ref")
    code2id <- .attrib_type_codes_to_ids_from_Ensembl_ftp(core_url, codes)
    if (toplevel) {
        keep_idx <- which(attrib_type_id == code2id[["toplevel"]])
        ids <- all_attribs[keep_idx, "seq_region_id"]
        toplevel <- seq_region_ids %in% ids
    } else {
        toplevel <- NULL
    }
    if (non_ref) {
        keep_idx <- which(attrib_type_id == code2id[["non_ref"]])
        ids <- all_attribs[keep_idx, "seq_region_id"]
        non_ref <- seq_region_ids %in% ids
    } else {
        non_ref <- NULL
    }
    list(toplevel=toplevel, non_ref=non_ref)
}

### Typical use:
###   core_url <- get_url_to_Ensembl_ftp_mysql_core("hsapiens_gene_ensembl")
###   coord_system_names <- c("chromosome", "scaffold", "lrg")
###   fetch_seq_regions_from_Ensembl_ftp(core_url, coord_system_names)
fetch_seq_regions_from_Ensembl_ftp <- function(core_url,
                                               coord_system_names=NULL,
                                               add.toplevel.col=FALSE,
                                               add.non_ref.col=FALSE)
{
    stopifnot(isTRUEorFALSE(add.toplevel.col),
              isTRUEorFALSE(add.non_ref.col))

    coord_systems <- fetch_default_coord_systems_from_Ensembl_ftp(core_url)
    Rtable <- "coord_system"
    colnames(coord_systems) <- paste(Rtable, colnames(coord_systems), sep=".")

    ## Even though we expect that column "coord_system.name" will always
    ## have unique values within the remaining rows, we want to be informed
    ## if this is not the case.
    coord_system.name <- coord_systems[ , "coord_system.name"]
    if (anyDuplicated(coord_system.name))
        warning(wmsg("column \"coord_system.name\" contains duplicates ",
                     "within the rows that have the default_version attrib"))

    if (!is.null(coord_system_names)) {
        ## Drop rows for which column "coord_system.name" is not
        ## in 'coord_system_names'.
        if (!is.character(coord_system_names))
            stop(wmsg("'coord_system_names' must be ",
                      "a character vector or NULL"))
        stop_if_not_primary_key(coord_system_names, "'coord_system_names'")
        keep_idx <- which(coord_system.name %in% coord_system_names)
        coord_systems <- S4Vectors:::extract_data_frame_rows(coord_systems,
                                                             keep_idx)
    }

    ## Fetch "seq_region" table and INNER JOIN with "coord_system" table.
    seq_regions <- fetch_table_from_Ensembl_ftp(core_url, "seq_region")
    Lcolname <- "coord_system_id"
    Rcolname <- paste(Rtable, Lcolname, sep=".")
    ans <- join_dfs(seq_regions, coord_systems, Lcolname, Rcolname)
    seq_region_ids <- ans[ , "seq_region_id"]

    ans$synonyms <- .fetch_synonyms_from_Ensembl_ftp(core_url, seq_region_ids)

    if (add.toplevel.col || add.non_ref.col) {
        cols <- .fetch_attribs_from_Ensembl_ftp(core_url,
                                                seq_region_ids,
                                                toplevel=add.toplevel.col,
                                                non_ref=add.non_ref.col)
        ans <- cbind(ans, S4Vectors:::delete_NULLs(cols),
                     stringsAsFactors=FALSE)
    }

    ans
}

