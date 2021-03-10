### =========================================================================
### Some low-level utilities to fetch data from Ensembl
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Obtain URL to the core db of a given species on the Ensembl FTP servers
###

.ENSEMBL_FTP_PUB_URL <- "ftp://ftp.ensembl.org/pub/"
.ENSEMBL_FTP_PUB_GRCH37_URL <- "ftp://ftp.ensembl.org/pub/grch37/"
.ENSEMBLGENOMES_FTP_PUB_URL <- "ftp://ftp.ensemblgenomes.org/pub/"
.ENSEMBL_FTP_RELEASE_PREFIX <- "release-"

check_species_info <- function(species_info)
{
    stopifnot(is.list(species_info))
    fields <- names(species_info)
    stopifnot(!is.null(fields), is_primary_key(fields),
              all(c("species", "Ensembl_release") %in% fields))
}

### 'division' must be NA or one of the Ensembl Genomes divisions i.e.
### "bacteria", "fungi", "metazoa", "plants", or "protists".
.get_Ensembl_FTP_top_url <- function(division=NA, use.grch37=FALSE)
{
    if (!is_single_value(division))
        stop(wmsg("'division' must be a single value"))
    if (!isTRUEorFALSE(use.grch37))
        stop(wmsg("'use.grch37' must be TRUE or FALSE"))
    if (!is.na(division)) {
        if (!is.character(division) || division == "")
            stop(wmsg("'division' must be a single non-empty string or NA"))
        if (use.grch37)
            stop(wmsg("'division' and 'use.grch37' cannot both be specified"))
        top_url <- paste0(.ENSEMBLGENOMES_FTP_PUB_URL, division, "/")
    } else if (use.grch37) {
        top_url <- .ENSEMBL_FTP_PUB_GRCH37_URL
    } else {
        top_url <- .ENSEMBL_FTP_PUB_URL
    }
    top_url
}

### The keys are FTP URLs to Ensembl division top-level directories e.g.
###   "ftp://ftp.ensembl.org/pub/"
###   "ftp://ftp.ensembl.org/pub/grch37/"
###   "ftp://ftp.ensemblgenomes.org/pub/plants/"
### etc...
.Ensembl_FTP_cached_releases <- new.env(parent=emptyenv())

list_Ensembl_FTP_releases <- function(division=NA, use.grch37=FALSE,
                                      as.subdirs=FALSE)
{
    if (!isTRUEorFALSE(as.subdirs))
        stop(wmsg("'as.subdirs' must be TRUE or FALSE"))
    top_url <- .get_Ensembl_FTP_top_url(division=division,
                                        use.grch37=use.grch37)
    releases <- .Ensembl_FTP_cached_releases[[top_url]]
    if (is.null(releases)) {
        top_files <- list_ftp_dir(top_url)
        nc <- nchar(.ENSEMBL_FTP_RELEASE_PREFIX)
        prefixes <- substr(top_files, 1L, nc)
        releases <- top_files[prefixes == .ENSEMBL_FTP_RELEASE_PREFIX]
        releases <- substr(releases, nc + 1L, nchar(releases))
        releases <- sort(as.integer(releases))
        .Ensembl_FTP_cached_releases[[top_url]] <- releases
    }
    if (as.subdirs)
        releases <- paste0(.ENSEMBL_FTP_RELEASE_PREFIX, releases)
    releases
}

### Returns a single string with the "species_info" attribute (named list)
### on it.
.predict_core_db_in_Ensembl_FTP_grch37 <- function(species=NA, release=NA)
{
    ## User-supplied 'species' will be ignored.
    warn <- TRUE
    if (missing(species))
        warn <- FALSE
    if (warn && is_single_value(species)) {
        if (is.na(species)) {
            ok <- TRUE
        } else if (is.character(species)) {
            ok <- species == "" ||
                  grepl("human", species, ignore.case=TRUE) ||
                  grepl("homo", species, ignore.case=TRUE) ||
                  grepl("sapiens", species, ignore.case=TRUE) ||
                  grepl("GRCh37", species, ignore.case=TRUE)
        } else {
            ok <- FALSE
        }
        warn <- !ok
    }
    if (warn)
        warning(wmsg("you've set 'use.grch37' to TRUE ",
                     "so 'species' was ignored"))
    if (is.na(release)) {
        release <- tail(list_Ensembl_FTP_releases(), n=1L)  # latest
    } else {
        release <- as.integer(release)
    }
    core_db <- paste0("homo_sapiens_core_", release, "_37")
    attr(core_db, "species_info") <- list(
        name="Human",
        species="homo_sapiens",
        division="EnsemblVertebrates",
        Ensembl_release=release,
        taxonomy_id=9606L,
        assembly="GRCh37",
        assembly_accession="GCF_000001405.13",
        core_db=core_db
    )
    core_db
}

### The species index file can be used with Ensembl releases >= 96 and
### Ensembl Genomes releases >= 22.
use_species_index_from_Ensembl_FTP <- function(release=NA, division=NA,
                                               use.grch37=FALSE)
{
    if (!is_single_value(release))
        stop(wmsg("'release' must be a single value"))
    if (!is_single_value(division))
        stop(wmsg("'division' must be a single value"))
    if (!isTRUEorFALSE(use.grch37))
        stop(wmsg("'use.grch37' must be TRUE or FALSE"))
    if (use.grch37) {
        if (is.na(division))
            return(FALSE)
        stop(wmsg("'division' and 'use.grch37' cannot both be specified"))
    }
    if (is.na(release))
        return(TRUE)
    release <- as.integer(release)
    if (is.na(division))
        return(release >= 96L)
    release >= 22L
}

.get_Ensembl_FTP_species_index_url <- function(release=NA, division=NA)
{
    if (!is_single_value(release))
        stop(wmsg("'release' must be a single value"))
    top_url <- .get_Ensembl_FTP_top_url(division=division)
    if (is.na(division)) {
        ## Available in Ensembl release 96 (March 2019) and above.
        species_file <- "species_EnsemblVertebrates.txt"
    } else {
        ## Available in Ensembl Genomes release 17 (Feb 2013) and above.
        ## However the current format is used only since release 22 (March
        ## 2014). See .load_or_fetch_species_index_from_url() below.
        species_file <- switch(division,
            bacteria="species_EnsemblBacteria.txt",
            fungi="species_EnsemblFungi.txt",
            metazoa="species_EnsemblMetazoa.txt",
            plants="species_EnsemblPlants.txt",
            protists="species_EnsemblProtists.txt",
            stop("Invalid division: ", division, "\n  ",
                 wmsg("Must be one of NA (stands for the main Ensembl ",
                      "division), \"bacteria\", \"fungi\", \"metazoa\", ",
                      "\"plants\", or \"protists\"."))
        )
    }
    if (!is.na(release)) {
        subdir <- paste0(.ENSEMBL_FTP_RELEASE_PREFIX, release)
    } else if (is.na(division)) {
        ## Get latest release.
        subdir <- tail(list_Ensembl_FTP_releases(as.subdirs=TRUE), n=1L)
    } else {
        subdir <- "current"
    }
    paste0(top_url, subdir, "/", species_file)
}

.fetch_species_index_from_url <- function(url)
{
    species_index <- fetch_table_from_url(url, header=TRUE)
    species_index_ncol <- ncol(species_index)
    if (species_index_ncol != 15L && species_index_ncol != 16L)
        stop(wmsg(url, " does not contain the expected fields"))

    ## This is the format used in Ensembl releases >= 96 and Ensembl
    ## Genomes releases >= 22.
    expected_colnames <- c(
        "name", "species", "division", "taxonomy_id",
        "assembly", "assembly_accession", "genebuild", "variation",
        "pan_compara", "peptide_compara", "genome_alignments",
        "other_alignments", "core_db", "species_id"
    )
    ## The "microarray" field was added in Ensembl 103 and Ensembl Genomes 50.
    if (species_index_ncol == 16L)
        expected_colnames <- append(expected_colnames, "microarray", after=8L)

    ## Note that Ensembl species index files are broken: the header line
    ## specifies 14 fields separated by 13 tabs BUT each line of data
    ## contains 14 tabs! The last tab is an additional tab placed
    ## at the end of the line i.e. after the 14th value in the line.
    ## For read.table() this means that there are actually 15 values
    ## and that the last value is missing! As a consequence the colnames
    ## on the returned data frame are completely messed up!
    expected_messed_up_colnames <- c("row.names", expected_colnames)
    expected_messed_up_colnames[2L] <- "X.name"
    if (!identical(colnames(species_index), expected_messed_up_colnames))
        stop(wmsg(url, " does not contain the expected fields"))
    colnames(species_index) <- c(expected_colnames,
                                 paste0("V", species_index_ncol))

    ## The last column should be filled with NAs. Drop it.
    stopifnot(all(is.na(species_index[[species_index_ncol]])))
    species_index[-species_index_ncol]
}

.cached_species_index <- new.env(parent=emptyenv())

.load_or_fetch_species_index_from_url <- function(url)
{
    species_index <- .cached_species_index[[url]]
    if (is.null(species_index)) {
        species_index <- .fetch_species_index_from_url(url)
        .cached_species_index[[url]] <- species_index
    }
    species_index
}

fetch_species_index_from_Ensembl_FTP <- function(release=NA, division=NA)
{
    url <- .get_Ensembl_FTP_species_index_url(release=release,
                                              division=division)
    .load_or_fetch_species_index_from_url(url)
}

.stop_on_ambiguous_lookup <- function(species_index, idx, max_print,
                                      species, column, url)
{
    show_matching_entries <- function(species_index, idx) {
        cat("\n  Matching entries")
        truncate <- length(idx) > max_print
        if (truncate) {
            cat(" (showing the first", max_print, "only)")
            idx <- head(idx, n=max_print)
        }
        cat(":\n")
        ## We drop columns that are totally uninteresting/meaningless in
        ## the context of species lookup.
        drop_columns <- c(
            "variation", "pan_compara", "peptide_compara",
            "genome_alignments", "other_alignments", "species_id"
        )
        m <- as.matrix(drop_cols(species_index[idx, ], drop_columns))
        if (truncate) {
            ellipsis <- "..."
            m <- rbind(m, rep.int(ellipsis, ncol(m)))
            rownames(m)[max_print + 1L] <- ellipsis
        }
        rownames(m) <- paste0("    ", rownames(m))
        print(m, quote=FALSE, right=TRUE)
    }
    on.exit(show_matching_entries(species_index, idx))
    stop("\n  ", wmsg("Found ", length(idx), " matches (case insensitive) ",
              "for \"", species, "\" in \"", column, "\" column of ",
              "Ensembl species index file:"),
              "\n    ", url)
}

### Matches a 'species' supplied as the name of a BioMart dataset (e.g.
### "hsapiens" or "hsapiens_gene_ensembl") to its corresponding core db.
.lookup_abbrev_species_in_core_dbs <- function(species, core_dbs)
{
    trimmed_core_dbs <- sub("_core_.*$", "", core_dbs)
    abbrev_core_dbs <- sub("^(.)[^_]*_", "\\1", trimmed_core_dbs)
    if (species == "mfuro_gene_ensembl") {
        abbrev_species <- "mputorius_furo"
    } else {
        abbrev_species <- strsplit(species, "_", fixed=TRUE)[[1L]][[1L]]
    }
    which(abbrev_core_dbs == abbrev_species)
}

.lookup_species <- function(species, species_index, url)
{
    .normalize_species <- function(species) chartr(" ", "_", tolower(species))

    ## Exact match (case insensitive).
    search_columns <- c("name", "species", "taxonomy_id",
                        "assembly", "assembly_accession", "core_db")
    for (column in search_columns) {
        if (column %in% c("name", "taxonomy_id")) {
            target <- species
        } else {
            target <- .normalize_species(species)
        }
        idx <- which(tolower(species_index[ , column]) %in% target)
        if (length(idx) == 1L)
            return(idx)
        if (length(idx) > 1L)
            .stop_on_ambiguous_lookup(species_index, idx, 4L,
                                      species, column, url)
    }

    ## Matches a 'species' supplied as the name of a BioMart dataset (e.g.
    ## "hsapiens_gene_ensembl" or "hsapiens") to the "core_db" column.
    idx <- .lookup_abbrev_species_in_core_dbs(.normalize_species(species),
                                              species_index[ , "core_db"])
    if (length(idx) == 1L)
        return(idx)
    if (length(idx) > 1L)
        .stop_on_ambiguous_lookup(species_index, idx, 4L,
                                  species, "core_db", url)

    ## Fuzzy match (grep()-based).
    search_columns <- c("name", "species", "assembly", "core_db")
    for (column in search_columns) {
        if (column %in% "name") {
            target <- species
        } else {
            target <- .normalize_species(species)
        }
        idx <- grep(target, species_index[ , column], ignore.case=TRUE)
        if (length(idx) == 1L)
            return(idx)
        if (length(idx) > 1L)
            .stop_on_ambiguous_lookup(species_index, idx, 4L,
                                      species, column, url)
    }

    ## We tried hard and failed!
    stop(wmsg("Couldn't find \"", species, "\" in ",
              "Ensembl species index file:"),
         "\n    ", url)
}

### Returns a single string with the "species_info" attribute (named list)
### on it.
.find_core_db_in_Ensembl_FTP_species_index <- function(species,
                                                       release=NA, division=NA)
{
    stopifnot(isSingleString(species))
    url <- .get_Ensembl_FTP_species_index_url(release=release,
                                              division=division)
    species_index <- .load_or_fetch_species_index_from_url(url)
    idx <- .lookup_species(species, species_index, url)
    core_db <- species_index[idx, "core_db"]
    if (is.na(release)) {
        release <- tail(list_Ensembl_FTP_releases(), n=1L)  # latest
    } else {
        release <- as.integer(release)
    }
    attr(core_db, "species_info") <- list(
        name=species_index[idx, "name"],
        species=species_index[idx, "species"],
        division=species_index[idx, "division"],
        Ensembl_release=release,
        taxonomy_id=species_index[idx, "taxonomy_id"],
        assembly=species_index[idx, "assembly"],
        assembly_accession=species_index[idx, "assembly_accession"],
        core_db=core_db
    )
    core_db
}

### 'division' must be NA or one of the Ensembl Genomes divisions i.e.
### "bacteria", "fungi", "metazoa", "plants", or "protists".
get_Ensembl_FTP_mysql_url <- function(release=NA, division=NA,
                                      use.grch37=FALSE)
{
    if (!is_single_value(release))
        stop(wmsg("'release' must be a single value"))
    top_url <- .get_Ensembl_FTP_top_url(division=division,
                                        use.grch37=use.grch37)
    if (!is.na(release)) {
        mysql_subdir <- paste0(.ENSEMBL_FTP_RELEASE_PREFIX, release, "/mysql")
    } else if (is.na(division) && !use.grch37) {
        mysql_subdir <- "current_mysql"
    } else {
        mysql_subdir <- "current/mysql"
    }
    paste0(top_url, mysql_subdir, "/")
}

### 'division' must be NA or one of the Ensembl Genomes divisions i.e.
### "bacteria", "fungi", "metazoa", "plants", or "protists".
get_Ensembl_FTP_gtf_url <- function(release=NA, division=NA,
                                    use.grch37=FALSE)
{
    if (!is_single_value(release))
        stop(wmsg("'release' must be a single value"))
    top_url <- .get_Ensembl_FTP_top_url(division=division,
                                        use.grch37=use.grch37)
    if (!is.na(release)) {
        gtf_subdir <- paste0(.ENSEMBL_FTP_RELEASE_PREFIX, release, "/gtf")
    } else if (is.na(division) && !use.grch37) {
        gtf_subdir <- "current_gtf"
    } else {
        gtf_subdir <- "current/gtf"
    }
    paste0(top_url, gtf_subdir, "/")
}

### The keys are FTP URLs to "mysql" directories e.g.
###   "ftp://ftp.ensembl.org/pub/current_mysql/"
###   "ftp://ftp.ensembl.org/pub/release-98/mysql/"
###   "ftp://ftp.ensemblgenomes.org/pub/bacteria/current/mysql/"
###   "ftp://ftp.ensemblgenomes.org/pub/plants/release-45/mysql/"
### etc...
.Ensembl_FTP_cached_core_dbs <- new.env(parent=emptyenv())

.list_Ensembl_FTP_core_dbs <- function(mysql_url, release=NA)
{
    stopifnot(isSingleString(mysql_url))

    pattern <- "_core_"
    core_dbs <- .Ensembl_FTP_cached_core_dbs[[mysql_url]]
    if (is.null(core_dbs)) {
        subdirs <- list_ftp_dir(mysql_url, subdirs.only=TRUE)
        core_dbs <- subdirs[grep(pattern, subdirs, fixed=TRUE)]
        .Ensembl_FTP_cached_core_dbs[[mysql_url]] <- core_dbs
    }

    if (!is.na(release))
        pattern <- paste0(pattern, release, "_")
    core_dbs[grep(pattern, core_dbs, fixed=TRUE)]
}

### Returns a single string with the "species_info" attribute (named list)
### on it.
.find_core_db_in_Ensembl_FTP_core_dbs <-
    function(mysql_url, species, release=NA)
{
    stopifnot(isSingleString(species))
    core_dbs <- .list_Ensembl_FTP_core_dbs(mysql_url, release=release)
    species2 <- chartr(" ", "_", tolower(species))
    idx <- .lookup_abbrev_species_in_core_dbs(species2, core_dbs)
    if (length(idx) != 1L)
        stop(wmsg("found 0 or more than 1 subdir for \"", species,
                  "\" species at ", mysql_url))
    core_db <- core_dbs[[idx]]
    if (is.na(release)) {
        release <- tail(list_Ensembl_FTP_releases(), n=1L)  # latest
    } else {
        release <- as.integer(release)
    }
    attr(core_db, "species_info") <- list(
        species=sub("_core_.*$", "", core_db),
        Ensembl_release=release,
        core_db=core_db
    )
    core_db
}

### 'division' must be NA or one of the Ensembl Genomes divisions i.e.
### "bacteria", "fungi", "metazoa", "plants", or "protists".
### Return URL to Ensemble Core DB (FTP access) in a single string with
### the "species_info" attribute (named list) on it.
get_Ensembl_FTP_core_db_url <- function(species, release=NA, division=NA,
                                        use.grch37=FALSE)
{
    mysql_url <- get_Ensembl_FTP_mysql_url(release=release,
                                           division=division,
                                           use.grch37=use.grch37)
    if (use.grch37) {
        core_db <- .predict_core_db_in_Ensembl_FTP_grch37(
                                           species,
                                           release=release)
    } else {
        if (!isSingleString(species) || species == "")
            stop(wmsg("'species' must be a single non-empty string"))
        use_species_index <- use_species_index_from_Ensembl_FTP(
                                           release=release,
                                           division=division,
                                           use.grch37=use.grch37)
        if (use_species_index) {
            ## The new way!
            core_db <- .find_core_db_in_Ensembl_FTP_species_index(
                                           species,
                                           release=release,
                                           division=division)
        } else {
            ## The old way. Dumb, unreliable, and slow!
            ## Only works if 'species' is supplied as the name of a BioMart
            ## dataset e.g. "celegans_gene_ensembl" (or "celegans").
            ## Ridiculously stiff and outdated!
            core_db <- .find_core_db_in_Ensembl_FTP_core_dbs(
                                           mysql_url,
                                           species, release=release)
        }
    }
    species_info <- attr(core_db, "species_info")
    check_species_info(species_info)
    core_db_url <- paste0(mysql_url, core_db, "/")
    attr(core_db_url, "species_info") <- species_info
    core_db_url
}


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
### Fetch core db data from the Ensembl FTP server (via RCurl)
###
### The utililities in this section use RCurl or just utils::download.file()
### for getting core db stuff directly from the Ensembl FTP server. They can
### access stuff that is not available thru biomaRt like for example the
### lengths of the sequences in the reference genome associated with a
### particular dataset and Ensembl release (e.g. for dataset
### "hsapiens_gene_ensembl" in release "64").
###
### Note that querying the Ensembl MySQL server (via RMariaDB) would
### probably be a better way to do this.
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

### Uses the utils::download.file() + utils::read.table() method.
fetch_table_from_Ensembl_FTP <- function(core_db_url, table,
                                         full.colnames=FALSE, nrows=-1L)
{
    columns <- .ENSEMBLDB_COLUMNS[[table]]
    if (is.null(columns)) {
        warning(wmsg("unknown table: ", table, " (download might fail or ",
                     "the returned data frame will have automatic colnames)"))
    } else if (full.colnames) {
        columns <- paste(table, columns, sep=".")
    }
    url <- paste0(core_db_url, table, ".txt.gz")
    ans <- fetch_table_from_url(url, colnames=columns, nrows=nrows)
    if (full.colnames && is.null(columns))
        colnames(ans) <- paste(table, colnames(ans), sep=".")
    ans
}

fetch_default_coord_systems_from_Ensembl_FTP <- function(core_db_url)
{
    coord_system <- fetch_table_from_Ensembl_FTP(core_db_url, "coord_system")

    ## Drop rows that do not have the default_version attrib.
    keep_idx <- grep("default_version", coord_system[ , "attrib"], fixed=TRUE)
    ans <- S4Vectors:::extract_data_frame_rows(coord_system, keep_idx)

    ## Even though we expect that column "name" will always have unique
    ## values within the remaining rows, we want to be informed if this
    ## is not the case.
    coord_system.name <- ans[ , "name"]
    if (anyDuplicated(coord_system.name))
        warning(wmsg("column \"coord_system.name\" contains duplicates ",
                     "within the rows that have the default_version attrib"))
    ans
}

### Retrieves all synonyms for the supplied sequence ids.
### This is done via tables "seq_region_synonym" and "external_db".
### Return a list parallel to 'seq_region_ids' where each list element is a
### named character vector containing all the synonyms for the corresponding
### sequence id. The names on the character vector indicate the origin of
### the synonym i.e. the external db where it's used (e.g. "INSDC",
### "RefSeq_genomic", "UCSC", "ensembl_internal_synonym", etc..)
.fetch_synonyms_from_Ensembl_FTP <- function(core_db_url, seq_region_ids)
{
    stopifnot(is_primary_key(seq_region_ids))

    all_synonyms <- fetch_table_from_Ensembl_FTP(core_db_url,
                                                 "seq_region_synonym")
    keep_idx <- which(all_synonyms[ , "seq_region_id"] %in% seq_region_ids)
    synonyms <- S4Vectors:::extract_data_frame_rows(all_synonyms, keep_idx)

    external_dbs <- fetch_table_from_Ensembl_FTP(core_db_url, "external_db")
    synonyms <- join_dfs(synonyms, external_dbs,
                         "external_db_id", "external_db_id",
                         keep.Rcol=TRUE)  # we'll drop it below

    f <- factor(synonyms[ , "seq_region_id"], levels=seq_region_ids)
    unlisted_ans <- setNames(synonyms[ , "synonym"], synonyms[ , "db_name"])
    split(unlisted_ans, f)
}

.attrib_type_codes_to_ids_from_Ensembl_FTP <- function(core_db_url, codes)
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
    attrib_type <- fetch_table_from_Ensembl_FTP(core_db_url, "attrib_type",
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

.external_db_names_to_ids_from_Ensembl_FTP <- function(core_db_url, db_names)
{
    if (!is.character(db_names))
        stop(wmsg("'db_names' must be a character vector"))
    if (anyDuplicated(db_names))
        stop(wmsg("'db_names' cannot contain duplicates"))
    external_db <- fetch_table_from_Ensembl_FTP(core_db_url, "external_db")
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
.fetch_attribs_from_Ensembl_FTP <- function(core_db_url, seq_region_ids,
                                            toplevel=FALSE, non_ref=FALSE)
{
    all_attribs <- fetch_table_from_Ensembl_FTP(core_db_url,
                                                "seq_region_attrib")
    attrib_type_id <- all_attribs[ , "attrib_type_id"]
    codes <- c("toplevel", "non_ref")
    code2id <- .attrib_type_codes_to_ids_from_Ensembl_FTP(core_db_url, codes)
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

### This is the workhorse behind getChromInfoFromEnsembl().
### Typical use:
###   core_db_url <- get_Ensembl_FTP_core_db_url("hsapiens")
###   fetch_seq_regions_from_Ensembl_FTP(core_db_url)
fetch_seq_regions_from_Ensembl_FTP <- function(core_db_url,
                                               add.toplevel.col=FALSE,
                                               add.non_ref.col=FALSE)
{
    stopifnot(isTRUEorFALSE(add.toplevel.col),
              isTRUEorFALSE(add.non_ref.col))

    coord_systems <- fetch_default_coord_systems_from_Ensembl_FTP(core_db_url)

    ## Fetch table "seq_region".
    seq_regions <- fetch_table_from_Ensembl_FTP(core_db_url, "seq_region")

    ## INNER JOIN table "seq_region" with table "coord_system".
    Rtable <- "coord_system"
    colnames(coord_systems) <- paste(Rtable, colnames(coord_systems), sep=".")
    Lcolumn <- "coord_system_id"
    Rcolumn <- paste(Rtable, Lcolumn, sep=".")
    ans <- join_dfs(seq_regions, coord_systems, Lcolumn, Rcolumn)
    seq_region_ids <- ans[ , "seq_region_id"]

    ans$synonyms <- .fetch_synonyms_from_Ensembl_FTP(core_db_url,
                                                     seq_region_ids)

    if (add.toplevel.col || add.non_ref.col) {
        cols <- .fetch_attribs_from_Ensembl_FTP(core_db_url,
                                                seq_region_ids,
                                                toplevel=add.toplevel.col,
                                                non_ref=add.non_ref.col)
        ans <- cbind(ans, S4Vectors:::delete_NULLs(cols),
                     stringsAsFactors=FALSE)
    }

    ans
}

