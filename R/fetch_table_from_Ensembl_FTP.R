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
### Note that fetching is done with fetch_table_from_url() which uses a
### 2-step approach (utils::download.file() && utils::read.table()).
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

### 'core_db_url' must be the full URL to a core DB directory located
### on the Ensembl FTP server e.g.
### "ftp://ftp.ensembl.org/pub/release-99/mysql/mus_musculus_core_99_38/"
### The "ftp://" part and trailing slash are mandatory!
### Use get_Ensembl_FTP_core_db_url() defined in Ensembl-utils.R to obtain
### such URL for a given species/release/division programmatically.
fetch_table_from_Ensembl_FTP <- function(core_db_url, table,
                                         full.colnames=FALSE, nrows=-1L)
{
    columns <- .ENSEMBLDB_COLUMNS[[table]]
    if (is.null(columns)) {
        warning(wmsg("unknown table: ", table, " (download might fail or ",
                     "the returned data frame will have automatic colnames)"),
                immediate.=TRUE)
    } else if (full.colnames) {
        columns <- paste(table, columns, sep=".")
    }
    url <- paste0(core_db_url, table, ".txt.gz")
    ans <- fetch_table_from_url(url, colnames=columns, nrows=nrows)
    if (full.colnames && is.null(columns))
        colnames(ans) <- paste(table, colnames(ans), sep=".")
    ans
}

