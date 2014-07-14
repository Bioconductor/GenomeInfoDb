### =========================================================================
### Convenience wrappers to the seqlevels() getter and setter
### -------------------------------------------------------------------------

keepSeqlevels <- function(x, value, ...)
{
    value <- unname(value)
    if (any(nomatch <- !value %in% seqlevels(x)))
        warning("invalid seqlevels '", 
                paste(value[nomatch], collapse=","), "' were ignored")
    if (is(x, "BSgenome"))
        stop("seqlevels cannot be dropped from a BSgenome object")
    force <- is(x, "Seqinfo")
    seqlevels(x, force=!force) <- value[!nomatch]
    x
}

dropSeqlevels <- function(x, value, ...)
{
    value <- unname(value)
    if (any(nomatch <- !value %in% seqlevels(x)))
        warning("invalid seqlevels '", 
                paste(value[nomatch], collapse=","), "' were ignored")
    if (is(x, "BSgenome"))
        stop("seqlevels cannot be dropped from a BSgenome object")
    force <- is(x, "Seqinfo")
    seqlevels(x, force=!force) <- seqlevels(x)[!seqlevels(x) %in% value]
    x
}

renameSeqlevels <- function(x, value, ...)
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
            idx <- match(nms, level)
            idx <- idx[complete.cases(idx)]
            value <- replace(level, idx, value[seq_along(idx)])
        } 
    } 
    seqlevels(x) <- value  
    x 
}

## Currently applies to TranscriptDb only.
restoreSeqlevels <- function(x, ...)
{
    seqlevels0(x) 
    x
}

keepStandardChromosomes <- function(x, species=NULL)
{
    style <- seqlevelsStyle(x)
    
    ori_seqlevels <- seqlevels(x)
    
    if(missing(species))
    {
        ## try to guess species here
        ## 1- get names of all species supported by GenomeInfodb
        ## 2- extract seqlevels for all species for current style
        ## 3- intersect, find how many species have some seqlevels as original?
        ## 4- if max(length(intersect)) is found once - then predict species
        allspecies <- names(genomeStyles())
    
        allseqlevels <- lapply(allspecies, function(x) 
            tryCatch(extractSeqlevels(x, style), error=function(e) NULL))
                
        mres <- sapply(allseqlevels, function(x) 
            length(intersect(x, ori_seqlevels)))
        
        if(length(which(mres == max(mres))) > 1)
            stop("Cannot determine standard chromosomes, Specify species arg")
        
        species <- allspecies[which(mres==max(mres))]
    } 
    
    standard_chromosomes <- extractSeqlevels(species, style)
    x <- keepSeqlevels(x,standard_chromosomes)
    return(x)
    
}

