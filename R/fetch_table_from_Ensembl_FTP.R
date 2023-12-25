### =========================================================================
### fetch_table_from_Ensembl_FTP()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ensembl Core schema (only for the tables we care about)
###
### List of Ensembl public MySQL servers / ports
###   https://www.ensembl.org/info/data/mysql.html
### Ensembl Core Schema:
###   https://www.ensembl.org/info/docs/api/core/core_schema.html

### Fundamental Tables

### 'attrib_type' table.
### (331 rows in homo_sapiens_core_99_38)
.ENSEMBLDB_ATTRIB_TYPE_COLUMNS <- c(
    "attrib_type_id",         # primary key
    "code",
    "name",
    "description"
)

### Assembly Tables

### 'seq_region' table.
### (268443 rows in homo_sapiens_core_99_38)
.ENSEMBLDB_SEQ_REGION_COLUMNS <- c(
    "seq_region_id",          # primary key
    "name",
    "coord_system_id",        # ==> coord_system.coord_system_id
    "length"
)

### 'coord_system' table.
### (9 rows in homo_sapiens_core_99_38)
.ENSEMBLDB_COORD_SYSTEM_COLUMNS <- c(
    "coord_system_id",        # primary key
    "species_id",
    "name",
    "version",
    "rank",
    "attrib"
)

### 'seq_region_attrib' table.
### (5927 rows in homo_sapiens_core_99_38)
.ENSEMBLDB_SEQ_REGION_ATTRIB_COLUMNS <- c(
    "seq_region_id",          # ==> seq_region.seq_region_id
    "attrib_type_id",         # ==> attrib_type.attrib_type_id
    "value"
)

### 'seq_region_synonym' table.
### (2360 rows in homo_sapiens_core_99_38)
.ENSEMBLDB_SEQ_REGION_SYNONYM_COLUMNS <- c(
    "seq_region_synonym_id",  # primary key
    "seq_region_id",          # ==> seq_region.seq_region_id
    "synonym",
    "external_db_id"          # ==> external_db.external_db_id
)

### External References

### 'external_db' table.
### (446 rows in homo_sapiens_core_99_38)
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
### fetch_table_from_Ensembl_FTP()
###
### Fetch a table (.txt.gz file) from one of the core DB dumps located on
### the Ensembl FTP server (ftp.ensembl.org).
### This provides access to stuff that is not available thru biomaRt like
### for example the lengths of the sequences in the reference genome
### associated with a particular dataset and Ensembl release (e.g. for
### dataset "hsapiens_gene_ensembl" in release "64").
###
### TODO: Querying the Ensembl MySQL server (via RMariaDB) would probably
### be a better way to do this.
### Update (Feb 10, 2020): Some preliminary testing indicates that using
### RMariaDB to fetch full tables is actually significantly slower.
### For example, to fetch table "seq_region" from the homo_sapiens_core_99_38
### core DB (the table has 268443 rows and 4 columns):
###   o The utils::download.file() + utils::read.table() approach takes about
###     7 sec to download and parse the MySQL dump located at:
###         ftp://ftp.ensembl.org/pub/current_mysql/homo_sapiens_core_99_38/
###   o The RMariaDB method takes about 22 sec to retrieve the table from the
###     MySQL server at ensembldb.ensembl.org.
###

.please_install_missing_CRAN_pkgs <- function(pkgs, reader0)
{
    fmt <- paste0("Couldn't load %s. The %s needed %s. Please install %s ",
                  "with 'install.packages(%s)' and try again.")
    if (length(pkgs) == 1L) {
        what1 <- paste0("package ", pkgs)
        what2 <- "package is"
        what3 <- "it"
        what4 <- paste0("\"", pkgs, "\"")
    } else {
        what1 <- paste0("the following packages: ", paste0(pkgs, collapse=", "))
        what2 <- "packages are"
        what3 <- "them"
        what4 <- paste0("c(", paste0("\"", pkgs, "\"",  collapse=", "), ")")
    }
    if (reader0 == "auto") {
        why <- paste0("to fetch data from Ensembl FTP server ",
                      "for Ensembl releases older than 99")
    } else {
        why <- paste0("when 'reader=\"", reader0, "\"'")
    }
    sprintf(fmt, what1, what2, why, what3, what4)
}

### 'core_db_url' must be the full URL to a core DB directory located
### on the Ensembl FTP server e.g.
### ftp://ftp.ensembl.org/pub/release-99/mysql/mus_musculus_core_99_38/ or
### ftp://ftp.ensembl.org/pub/grch37/release-87/mysql/homo_sapiens_core_87_37/
### The "ftp://" part and trailing slash are mandatory!
### Use get_Ensembl_FTP_core_db_url() defined in Ensembl-utils.R to obtain
### such URL for a given species/release/division programmatically.
fetch_table_from_Ensembl_FTP <-
    function(core_db_url, table, full.colnames=FALSE, nrows=-1L,
             reader=c("auto", "read.table", "fread"))
{
    reader0 <- reader <- match.arg(reader)
    columns <- .ENSEMBLDB_COLUMNS[[table]]
    if (!is.null(columns) && full.colnames)
        columns <- paste(table, columns, sep=".")
    url <- paste0(core_db_url, table, ".txt.gz")
    if (reader == "auto") {
        ## Table dumps from Ensembl release < 99 can contain inline \r
        ## characters (carriage returns) which break utils::read.table().
        ## However data.table::fread() seems to be able to handle them.
        ## See https://github.com/Bioconductor/GenomeInfoDb/issues/98.
        release <- sub("_.*$", "", sub("^.*_core_", "", core_db_url))
        release <- suppressWarnings(as.integer(release))
        if (is.na(release) || release >= 99) {
            reader <- "read.table"
        } else {
            reader <- "fread"
        }
    }
    if (reader == "read.table") {
        if (is.null(columns))
            warning(wmsg("unknown table: ", table, " (download might fail ",
                         "or the returned data frame will have automatic ",
                         "colnames)"),
                    immediate.=TRUE)
        ## fetch_table_from_url() downloads the full file before reading it.
        ans <- fetch_table_from_url(url, colnames=columns, nrows=nrows)
        if (is.null(columns) && full.colnames)
            colnames(ans) <- paste(table, colnames(ans), sep=".")
    } else {
        if (is.null(columns))
            stop(wmsg("unknown table: ", table))
        ## data.table::fread() needs R.utils to read compressed files.
        missing_pkgs <- character(0)
        if (!requireNamespace("R.utils", quietly=TRUE))
            missing_pkgs <- c(missing_pkgs, "R.utils")
        if (!requireNamespace("data.table", quietly=TRUE))
            missing_pkgs <- c(missing_pkgs, "data.table")
        if (length(missing_pkgs) != 0L) {
            errmsg <- .please_install_missing_CRAN_pkgs(missing_pkgs, reader0)
            stop(wmsg(errmsg))
        }
        ## data.table::fread() downloads the full file before reading it.
        ans <- data.table::fread(url, nrows=nrows, strip.white=FALSE,
                                      showProgress=FALSE)
        ans <- as.data.frame(ans)
        colnames(ans) <- columns
    }
    ans
}

