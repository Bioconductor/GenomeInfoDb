### =========================================================================
### Convenience wrappers to the seqlevels() getter and setter
### -------------------------------------------------------------------------

keepSeqlevels <- function(x, value,
                          pruning.mode=c("error", "coarse", "fine", "tidy"))
{
    if (is(x, "BSgenome"))
        stop("seqlevels cannot be dropped from a BSgenome object")
    if (!is.null(names(value)))
        warning("the names on the vector of supplied seqlevels were ignored")
    value <- as.character(value)
    is_invalid <- !(value %in% seqlevels(x))
    if (any(is_invalid))
        stop(wmsg("invalid seqlevels: ",
                  paste(unique(value[is_invalid]), collapse=", ")))
    seqlevels(x, pruning.mode=pruning.mode) <- value
    x
}

dropSeqlevels <- function(x, value,
                          pruning.mode=c("error", "coarse", "fine", "tidy"))
{
    if (is(x, "BSgenome"))
        stop("seqlevels cannot be dropped from a BSgenome object")
    if (!is.null(names(value)))
        warning("the names on the vector of supplied seqlevels were ignored")
    keep_seqlevels <- setdiff(seqlevels(x), as.character(value))
    seqlevels(x, pruning.mode=pruning.mode) <- keep_seqlevels
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

standardChromosomes <- function(x, species=NULL)
{
    ori_seqlevels <- seqlevels(x)
    if(!length(ori_seqlevels))
        return(character())

    ## guess at style
    guess <- .guessSpeciesStyle(ori_seqlevels)
    if (!any(is.na(guess)))
        style <- unique(guess$style)
    else
        return(character())
    if (length(style) > 1)
        style <- style[1]

    standard <- character()
    if (is.null(species)) {
        possible <- genomeStyles()
        ## extractSeqlevels will fail if no style-species match; 
        ## must check style match first
        standard <- unique(unlist(
            Map(function(name, data) {
                if (style %in% colnames(data))
                    intersect(ori_seqlevels, extractSeqlevels(name, style))
            }, name=names(possible), data=possible)))
        ## preserve order
        standard <- ori_seqlevels[ori_seqlevels %in% standard]
        if (!length(standard))
            stop(paste0("cannot determine standard chromosomes; ",
                        "try specifying 'species' argument"))
    } else {
        standard <- extractSeqlevels(species, style)
        standard <- intersect(ori_seqlevels,standard)
    }
    standard 
}

keepStandardChromosomes <-
    function(x, species=NULL, pruning.mode=c("error", "coarse", "fine", "tidy"))
{
    standard <- standardChromosomes(x, species)
    keepSeqlevels(x, standard, pruning.mode=pruning.mode)
}
