############################################################################
## Helpers for getting genus and species from the NCBI taxonomy ID
## And vice versa...
globalVariables(c("speciesMap"))

## modified function to use our own data.
.lookupSpeciesFromTaxId <- function(id){
    if (!exists("specData")) {
     data(specData, package = "GenomeInfoDb")
    }
    ## Then find matches
    g <- specData[,1] == id
    res <- specData[g,]
    if(dim(res)[1]<1) stop("Cannot find a species to match the requested taxonomy id. Please provide the genus and species manually.")
    if(dim(res)[1]==1) return(res[1,])
    if(dim(res)[1]>1){
        .tooLong <- function(x){
            splt <- unlist(strsplit(x,split=" "))
            if(length(splt) > 1){
                return(TRUE)
            }else{
                return(FALSE)
            }
        }
        tooLong = unlist(lapply(as.character(res$species), .tooLong))
        if(all(tooLong)){
            return(res[1,])
        }else{
            res <- res[!tooLong,]
            return(res[1,])
        }
    }
}

## modify the above to throw out species=NA and only keep the 1st result

## usage:
## .lookupSpeciesFromTaxId("10090")
## bigger test:
## res <- list()
## for(i in seq_along(taxIDs)){
##     message(paste0('now testing: ',taxIDs[i]))
##     res[[i]] <- .lookupSpeciesFromTaxId(taxIDs[i])
## }
 
## This indicates that I have about 1438 different organisms that I
## can run blast2GO for plus 38 organisms that I have NCBI GO data
## for.




## Gets a taxonomyId For a species name...
.taxonomyId <-
    function(species)
{
    if(is.na(species)){return(NA)}
    if (!exists("speciesMap"))
        data(speciesMap, package="GenomeInfoDb")
    species <- gsub(" {2,}", " ", species)
    species <- gsub(",", " ", species, fixed=TRUE)
    idx <- match(species, speciesMap$species)
    if (any(is.na(idx)))
        stop(sum(is.na(idx)), " unknown species: ",
             paste(sQuote(head(species[is.na(idx)])), collapse=" "))
    as.character(speciesMap$taxon[idx])
}

## usage:
## .taxonomyId("Homo sapiens")


## Checks to see if a tax Id is valid or not.
.checkForAValidTaxonomyId <- function(taxId){
## TODO: precompute the list of valid tax Ids
if (!exists("validTaxIds")) {
     data(validTaxIds, package = "GenomeInfoDb")
}
validTaxIds <- c(validTaxIds, NA_integer_)
if(!(taxId %in% validTaxIds)){
      stop(wmsg(paste0("The taxonomy Id you have provided (",taxId,")",
                       " is not in our list of valid Tax Ids.",
                       ". Please check to make sure that your tax ID",
                       " is really legitimate and if so, then please tell",
                       " us about it so that we can update our list."))) 
  }
}

## usage:
## .checkForAValidTaxonomyId(9606)




## NOTES:
## Right now the data sources for these different functions are different.  This is partly for performance reasons and partly for historical ones.  But today, we are just putting all this functionality into one place.
## Another thing that is different is that specData.rda has an origin that is slightly different than that of speciesMap.rda and validTaxIds.rda.  This is also historical, but all these resources came from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz (just at different times).  The code used to process specData.rrda was preserved (below) and can be modified to make all three resources as updates at a future time).




###################################################################3
## Code for processing the tax IDs directly from NCBI.
## So look here for a mapping file
## ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
## And grab out the names.dmp file
## and then do this to it to preprocess it:
.processTaxNamesFile <- function(){
    species <- read.delim('names.dmp',header = FALSE,sep = "|")
    colnames(species) <- c('tax_id','name_txt','unique_name','name_class')
    ## keep only 1st two cols
    species <- species[,c(1:2,4)]
    ## throw away tabs from second col
    species[[2]] <- gsub('\t','',species[[2]])
    ## And the third col
    species[[3]] <- gsub('\t','',species[[3]])
    ## throw away rows where the third column doesn't say 'scientific name'
    keep <- grepl('scientific name', species[[3]])
    species <- species[keep,1:2]
    
    ## split second column by first space:
    rawSpec <- species[[2]]
    spltSpec <- strsplit(rawSpec, split=" ")
    genusDat <- unlist(lapply(spltSpec, function(x){x[1]}))
    .getRest <- function(x){
        if(length(x) > 1){
            return(paste(x[2:length(x)], collapse=" "))
        }else{
            return(NA)
        }
    }
    speciesDat <- unlist(lapply(spltSpec, .getRest))
    specData <- data.frame(tax_id=species[[1]],
                           genus=genusDat,
                           species=speciesDat,
                           stringsAsFactors=FALSE)
    ## name columns
    save(specData, file='specData.rda')
}
## .processTaxNamesFile()



