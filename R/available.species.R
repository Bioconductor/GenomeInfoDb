### =========================================================================
### Helpers to map between genus, species and taxonomy ID 
### -------------------------------------------------------------------------

## In February 2017 the mapping files in GenomeInfoDb/data/ were moved to the 
## GenomeInfoDbData annotation package.

.lookupSpeciesFromTaxId <- function(id, all=FALSE) {
    if (!exists("specData")) {
     data(specData, package = "GenomeInfoDbData")
    }
    ## Then find matches
    g <- specData[,1] == id
    res <- specData[g,]
    if (dim(res)[1]<1) 
        stop(paste0("Cannot find a species to match the requested",
                    " taxonomy id. Please provide the genus and species",
                    " manually."))
    if (dim(res)[1] == 1) { 
        return(res[1,])
    } else if (dim(res)[1]>1) {
        .tooLong <- function(x){
            splt <- unlist(strsplit(x,split=" "))
            if(length(splt) > 1){
                return(TRUE)
            }else{
                return(FALSE)
            }
        }
        tooLong <- unlist(lapply(as.character(res$species), .tooLong))
        if (all(tooLong)) {
            return(res[1,])
        } else {
            res <- res[!tooLong,]
            if (all) {
                return(res)
            } else {
                return(res[1,])
            }
        }
    }
}

available.species <- function(){
    if (!exists("speciesMap"))
        data(speciesMap, package="GenomeInfoDbData")
    speciesMap
}

.getTaxonomyId <- function(species) {
    if (is.na(species)) {return(NA)} 
    if (!exists("speciesMap"))
        data(speciesMap, package="GenomeInfoDbData")
    species <- gsub(" {2,}", " ", species)
    species <- gsub(",", " ", species, fixed=TRUE)
    idx <- match(species, speciesMap$species)
    if (any(is.na(idx)))
        stop(sum(is.na(idx)), " unknown species: ",
             paste(sQuote(head(species[is.na(idx)])),
                   paste0("Please use 'available.species' to see viable",
                          " species names or tax Ids"),
                   collapse=" "))
    as.integer(speciesMap$taxon[idx])
}

.taxonomyId <- function(species){
    unlist(lapply(species, .getTaxonomyId))
}

.checkForAValidTaxonomyId <- function(taxId) {
    if (!exists("validTaxIds"))
        data(validTaxIds, package = "GenomeInfoDbData")
    validTaxIds <- c(validTaxIds, NA_integer_)
    if(!(taxId %in% validTaxIds)) {
          stop(wmsg(paste0("The taxonomy Id you have provided (",taxId,")",
                           " is not in our list of valid Tax Ids.",
                           " Please check to make sure that your tax ID",
                           " is really legitimate and if so, then please tell",
                           " us about it so that we can update our list."))) 
    }
}
