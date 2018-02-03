### =========================================================================
### Helpers to map between taxonomy ID, genus, and organism
### -------------------------------------------------------------------------

### In February 2017 the mapping files in GenomeInfoDb/data/ were moved to the
### GenomeInfoDbData annotation package.

.TAXONOMY_DB_cache <- new.env(parent=emptyenv())

### Return a data.frame with 3 columns: tax_id, genus, species.
### Number of rows: 1820543 (as of Jan 30, 2018).
### TODO: Rename specData dataset -> TAXONOMY_DB in GenomeInfoDbData.
loadTaxonomyDb <- function()
{
    ans <- try(get("TAXONOMY_DB", envir=.TAXONOMY_DB_cache, inherits=FALSE),
               silent=TRUE)
    if (!is(ans, "try-error"))
        return(ans)
    data(specData, package="GenomeInfoDbData", envir=.TAXONOMY_DB_cache)
    taxdb <- get("specData", envir=.TAXONOMY_DB_cache, inherits=FALSE)
    stopifnot(identical(colnames(taxdb), c("tax_id", "genus", "species")),
              is.integer(taxdb[["tax_id"]]),
              is.factor(taxdb[["genus"]]),
              is.character(taxdb[["species"]]))
    ## Replace NAs in the "species" column with emty strings.
    ## Shouldn't we clean the dataset in GenomeInfoDbData instead?
    taxdb[["species"]][is.na(taxdb[["species"]])] <- ""
    assign("TAXONOMY_DB", taxdb, envir=.TAXONOMY_DB_cache)
    taxdb
}

available.species <- function()
{
    .Deprecated("loadTaxonomyDb")
    loadTaxonomyDb()
}

### NOT exported but used in the GenomicFeatures package.
### Not vectorized.
lookup_organism_by_tax_id <- function(tax_id, all=FALSE)
{
    stopifnot(isSingleNumber(tax_id))
    taxdb <- loadTaxonomyDb()
    ## Find matches.
    idx <- which(taxdb[["tax_id"]] == tax_id)
    if (length(idx) == 0L)
        stop(wmsg("Cannot find an organism to match the requested ",
                  "taxonomy ID. Please provide the genus and species ",
                  "manually."))
    ans <- taxdb[idx, , drop=FALSE]
    if (nrow(ans) == 1L || all)
        return(ans)
    ## When nrow(ans) > 1 and 'all' is FALSE, we first reduce the set of
    ## entries to keep single word species only, then pick up the first
    ## entry.
    idx1 <- which(lengths(strsplit(ans[["species"]], split=" ")) == 1L)
    if (length(idx1) >= 1L)
        ans <- ans[idx1, , drop=FALSE]
    ans[1L, , drop=FALSE]
}

### NOT exported but used in the GenomicFeatures package.
### Not vectorized.
lookup_tax_id_by_organism <- function(organism)
{
    stopifnot(is.character(organism) || is.factor(organism),
              length(organism) == 1L)
    if (is.na(organism)) return(NA)
    taxdb <- loadTaxonomyDb()
    species <- taxdb[["species"]]
    organisms <- trimws(paste(taxdb[["genus"]], species))
    organism <- gsub(" {2,}", " ", organism)
    organism <- gsub(",", " ", organism, fixed=TRUE)
    idx <- match(organism, organisms)
    if (is.na(idx))
        stop(wmsg(organism, ": unknown organism. ",
                  "Please use 'loadTaxonomyDb()' to see viable ",
                  "genus/species and taxonomy IDs."))
    as.integer(taxdb[["tax_id"]][[idx]])
}

### NOT exported but used in the GenomicFeatures package.
### Vectorized.
check_tax_id <- function(tax_id)
{
    stopifnot(is.numeric(tax_id))
    if (!is.integer(tax_id))
        tax_id <- as.integer(tax_id)
    taxdb <- loadTaxonomyDb()
    bad_idx <- which(!(tax_id %in% taxdb[["tax_id"]]))
    if (length(bad_idx) != 0L) {
          bad_ids <- paste0(unique(tax_id[bad_idx]), collapse=", ")
          stop(wmsg("Unknown taxonomy IDs: ", bad_ids,
                    "\n\n  These taxonomy IDs are not in our list of valid ",
                    "taxonomy IDs. Please check to make sure that the ",
                    "supplied taxonomy IDs are legitimate and if so, then ",
                    "please tell us about it so that we can update our list."))
    }
}

