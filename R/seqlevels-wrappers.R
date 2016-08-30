### =========================================================================
### Convenience wrappers to the seqlevels() getter and setter
### -------------------------------------------------------------------------

keepSeqlevels <- function(x, value)
{
    value <- unname(value)
    if (any(nomatch <- !value %in% seqlevels(x)))
        warning("invalid seqlevels", 
                paste(sQuote(value[nomatch]), collapse=", "), "were ignored")
    if (is(x, "BSgenome"))
        stop("seqlevels cannot be dropped from a BSgenome object")
    force <- is(x, "Seqinfo")
    seqlevels(x, force=!force) <- value[!nomatch]
    x
}

dropSeqlevels <- function(x, value)
{
    value <- unname(value)
    if (any(nomatch <- !value %in% seqlevels(x)))
        warning("invalid seqlevels", 
                paste(sQuote(value[nomatch]), collapse=", "), "were ignored")
    if (is(x, "BSgenome"))
        stop("seqlevels cannot be dropped from a BSgenome object")
    force <- is(x, "Seqinfo")
    seqlevels(x, force=!force) <- seqlevels(x)[!seqlevels(x) %in% value]
    x
}

renameSeqlevels <- function(x, value)
{
    nms <- names(value)
    ## unnamed
    if (is.null(nms)) {
        if (length(value) != length(seqlevels(x)))
            stop("unnamed 'value' must be the same length as seqlevels(x)")
        names(value) <- seqlevels(x)
    ## named
    } else {
        if (any(nomatch <- !nms %in% seqlevels(x)))
            warning("invalid seqlevels ", 
                    paste(sQuote(nms[nomatch]), collapse=", "), " ignored")
        if (length(value) != length(seqlevels(x))) {
            level <- seqlevels(x)
            idx <- match(level, nms)
            level[!is.na(idx)] <- value[na.omit(idx)]
            value <- level
        } 
    } 
    seqlevels(x) <- value 
    x 
}

## Currently applies to TxDb only.
restoreSeqlevels <- function(x)
{
    seqlevels(x) <- seqlevels0(x) 
    x
}

keepStandardChromosomes <- function(x, species=NULL)
{
    ori_seqlevels <- seqlevels(x)
 
    if(length(ori_seqlevels)==0)
        return(x)
 
    ans <- .guessSpeciesStyle(ori_seqlevels)
 
    if(length(ans)==1){
        if(is.na(ans)){
            ## seqlevels did not match any organism's style - so
            ## drop all levels and return an empty object
            return(dropSeqlevels(x, seqlevels(x)))
        }
    }
 
    standard <- character(0)
    if(missing(species))
    {
        ## compatible species:
        mres <- Map(function(y, style) {
                    levels <- extractSeqlevels(y, style)
                    intersect(ori_seqlevels, levels)
                }, y = ans$species, style = ans$style) 

        ## seqlevels must be the same for all possible species, else fail
        standard <- unique(unlist(mres))
        if (length(unique(lengths(mres))) != 1L | 
            length(standard) != lengths(mres[1]) | 
            length(standard) == 0L)
            stop(paste0("cannot determine standard chromosomes;",
                 " try specifying 'species' arg"))
    } else {
        ## exact species:
        ## .guessSpeciesStyle() returned 'possible' species using weights
        idx <- which(ans$species == .normalize_organism(species))
        if (length(idx) == 0L)
            stop(paste0("'species' not found or not compatible with given ",
                 "seqlevels; see names(genomeStyles()) for valid species"))
        if (length(idx) > 1L)
            idx <- idx[1]
        standard <- extractSeqlevels(ans$species[idx], ans$style[idx])
        standard <- intersect(ori_seqlevels,standard)
    }

    seqlevels(x, force=TRUE) <- standard 
    x 
}

