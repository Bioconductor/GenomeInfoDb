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
            ## Interanally the seqlevels did not match any organism's style - so
            ## drop all levels and return an empty object
            return(dropSeqlevels(x, seqlevels(x)))
        }
    }
 
    style <- unique(ans$style)
    if (length(style) > 1L) {
        message(paste0("using the first of multiple styles matched: ",
                       paste(style, collapse=", ")))
        style <- style[1]
    }
    standard <- character(0)
 
    if(missing(species))
    {
        ## compatible species
        possible <- genomeStyles()
        idx <- sapply(possible, function(y) style %in% colnames(y))
        compatible <- names(possible)[idx]
        ## intersect seqlevels for compatible species with original
        mres <- sapply(compatible, function(y) {
            levels <- extractSeqlevels(y, style)
            intersect(levels, ori_seqlevels)
        }) 
 
        standard <- unique(unlist(mres))
        standard <- ori_seqlevels[which(ori_seqlevels %in% standard)]
        if(length(standard) == 0 )
            stop(paste0("cannot determine standard chromosomes;",
                 " try specifying 'species' arg"))
 
    }else{
        standard <- extractSeqlevels(species, style)
        standard <- intersect(ori_seqlevels,standard)
    }

    seqlevels(x, force=TRUE) <- standard 
    x 
}

