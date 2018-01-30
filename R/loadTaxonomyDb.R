### =========================================================================
### Helpers to map between taxonomy ID, genus, and species
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
lookup_organism_by_tax_id <- function(tax_id, all=FALSE)
{
    taxdb <- loadTaxonomyDb()
    ## Find matches.
    g <- taxdb[["tax_id"]] == tax_id
    ans <- taxdb[g, , drop=FALSE]
    if (nrow(ans) < 1) 
        stop(wmsg("Cannot find a species to match the requested ",
                  "taxonomy ID. Please provide the genus and species ",
                  "manually."))
    if (nrow(ans) == 1)
        return(ans)
    ## nrow(ans) > 1
    .tooLong <- function(x) {
        splt <- unlist(strsplit(x, split=" "))
        length(splt) > 1
    }
    tooLong <- unlist(lapply(ans[["species"]], .tooLong))
    if (all(tooLong))
        return(ans[1, ])
    ans <- ans[!tooLong, , drop=FALSE]
    if (!all)
        ans <- ans[1, ]
    ans
}

### NOT exported but used in the GenomicFeatures package.
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
check_tax_id <- function(tax_id)
{
    stopifnot(isSingleNumber(tax_id))
    if (!is.integer(tax_id))
        tax_id <- as.integer(tax_id)
    taxdb <- loadTaxonomyDb()
    if (!(tax_id %in% taxdb[["tax_id"]])) {
          stop(wmsg("The taxonomy ID you have provided (", tax_id, ") ",
                    "is not in our list of valid taxonomy IDs. ",
                    "Please check to make sure that your taxonomy ID ",
                    "is legitimate and if so, then please tell ",
                    "us about it so that we can update our list."))
    }
}

