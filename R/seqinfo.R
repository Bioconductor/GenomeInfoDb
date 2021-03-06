### =========================================================================
### The seqinfo() (and releated) generic getters and setters
### -------------------------------------------------------------------------


normarg_new2old <- function(new2old, new_N, old_N,
                            new_what="the supplied 'seqinfo'",
                            old_what="the current 'seqinfo'")
{
    if (!is.numeric(new2old))
        stop(wmsg("'new2old' must be NULL or an integer vector"))
    if (length(new2old) != new_N)
        stop(wmsg("when not NULL, 'new2old' must ",
                  "have the length of ", new_what))
    if (!is.integer(new2old))
        new2old <- as.integer(new2old)
    min_new2old <- suppressWarnings(min(new2old, na.rm=TRUE))
    if (min_new2old != Inf) {
        if (min_new2old < 1L)
            stop(wmsg("when not NULL, 'new2old' must ",
                      "contain positive values or NAs"))
        if (max(new2old, na.rm=TRUE) > old_N)
            stop(wmsg("'new2old' cannot contain values ",
                      "greater than the length of ", old_what))
    }
    if (any(duplicated(new2old) & !is.na(new2old)))
        stop(wmsg("non-NA values in 'new2old' must be unique"))
    new2old
}

### Validate and reverse the 'new2old' mapping.
.reverse_new2old <- function(new2old, new_N, old_N,
                             new_what="the supplied 'seqinfo'",
                             old_what="the current 'seqinfo'")
{
    if (is.null(new2old))
        return(NULL)
    new2old <- normarg_new2old(new2old, new_N, old_N, new_what, old_what)
    S4Vectors:::reverseIntegerInjection(new2old, old_N)
}

### The dangling seqlevels in 'x' are those seqlevels that the user wants to
### drop but are in use.
getDanglingSeqlevels <- function(x, new2old=NULL,
                            pruning.mode=c("error", "coarse", "fine", "tidy"),
                            new_seqlevels)
{
    pruning.mode <- match.arg(pruning.mode)
    if (!is.character(new_seqlevels) || any(is.na(new_seqlevels)))
        stop(wmsg("the supplied 'seqlevels' must be a character vector ",
                  "with no NAs"))
    if (is.null(new2old))
        return(character(0))
    new_N <- length(new_seqlevels)
    old_seqlevels <- seqlevels(x)
    old_N <- length(old_seqlevels)
    old2new <- .reverse_new2old(new2old, new_N, old_N,
                                new_what="the supplied 'seqlevels'",
                                old_what="the current 'seqlevels'")
    seqlevels_to_drop <- old_seqlevels[is.na(old2new)]
    seqlevels_in_use <- seqlevelsInUse(x)
    dangling_seqlevels <- intersect(seqlevels_to_drop, seqlevels_in_use)
    if (length(dangling_seqlevels) != 0L && pruning.mode == "error")
        stop(wmsg("The following seqlevels are to be dropped but are ",
                  "currently in use (i.e. have ranges on them): ",
                  paste(dangling_seqlevels, collapse = ", "), ".\n",
                  "Please use the 'pruning.mode' argument to control how ",
                  "to prune 'x', that is, how to remove the ranges in 'x' ",
                  "that are on these sequences. For example, do something ",
                  "like:"),
             "\n    seqlevels(x, pruning.mode=\"coarse\") <- new_seqlevels",
             "\n  or:",
             "\n    keepSeqlevels(x, new_seqlevels, pruning.mode=\"coarse\")",
             "\n  See ?seqinfo for a description of the pruning modes.")
    dangling_seqlevels
}

### Compute the new seqnames resulting from new seqlevels.
### Assumes that 'seqnames(x)' is a 'factor' Rle (which is true if 'x' is a
### GRanges or GAlignments object, but not if it's a GRangesList object),
### and returns a 'factor' Rle of the same length (and same runLength vector).
### Always used in the context of the seqinfo() setter i.e. 'new_seqlevels'
### comes from 'seqlevels(value)' where 'value' is the supplied Seqinfo object.
makeNewSeqnames <- function(x, new2old=NULL, new_seqlevels)
{
    ## Should never happen.
    stopifnot(is.character(new_seqlevels), all(!is.na(new_seqlevels)))
    new_N <- length(new_seqlevels)
    old_N <- length(seqlevels(x))
    x_seqnames <- seqnames(x)
    if (!is.null(new2old)) {
        old2new <- .reverse_new2old(new2old, new_N, old_N,
                                    new_what="the supplied Seqinfo object",
                                    old_what="seqinfo(x)")
        tmp <- runValue(x_seqnames)
        levels(tmp) <- new_seqlevels[old2new]
        runValue(x_seqnames) <- factor(as.character(tmp), levels=new_seqlevels)
        return(x_seqnames)
    }
    if (new_N >= old_N &&
        identical(new_seqlevels[seq_len(old_N)], seqlevels(x)))
    {
        levels(x_seqnames) <- new_seqlevels
        return(x_seqnames)
    }
    SEQLEVELS_ARE_NOT_THE_SAME <- c(
        "The seqlevels in the supplied Seqinfo object ",
        "are not the same as the seqlevels in 'x'. "
    )
    if (length(intersect(seqlevels(x), new_seqlevels)) == 0L)
        stop(wmsg(SEQLEVELS_ARE_NOT_THE_SAME,
                  "Please use the 'new2old' argument to specify the ",
                  "mapping between the formers and the latters."))
    if (setequal(seqlevels(x), new_seqlevels))
        stop(wmsg("The seqlevels in the supplied Seqinfo object ",
                  "are not in the same order as the seqlevels in 'x'. ",
                  "Please reorder the seqlevels in 'x' with:"),
             "\n\n",
             "    seqlevels(x) <- seqlevels(new_seqinfo)\n\n  ",
             wmsg("before calling the 'seqinfo()' setter."),
             "\n  ",
             wmsg("For any more complicated mapping between the new ",
                  "and old seqlevels (e.g. for a mapping that will ",
                  "result in the renaming of some seqlevels in 'x'), ",
                  "please use the 'new2old' argument."))
    if (all(seqlevels(x) %in% new_seqlevels))
        stop(wmsg(SEQLEVELS_ARE_NOT_THE_SAME,
                  "To map them to the seqlevels of the same name in 'x', ",
                  "the easiest way is to propagate them to 'x' with:"),
             "\n\n",
             "    seqlevels(x) <- seqlevels(new_seqinfo)\n\n  ",
             wmsg("before calling the 'seqinfo()' setter."),
             "\n  ",
             wmsg("For any more complicated mapping, please use ",
                  "the 'new2old' argument."))
    stop(wmsg(SEQLEVELS_ARE_NOT_THE_SAME,
              "To map them to the seqlevels of the same name in 'x', ",
              "the easiest way is to propagate them to 'x' with:"),
         "\n\n",
         "    seqlevels(x) <- seqlevels(new_seqinfo)\n\n  ",
         wmsg("before calling the 'seqinfo()' setter. ",
	      "Note that you might need to specify a pruning mode ",
              "(via the 'pruning.mode' argument) if this operation ",
              "will drop seqlevels that are in use in 'x'."),
         "\n  ",
         wmsg("For any more complicated mapping, please use ",
              "the 'new2old' argument."))
}

### Return -3L for "renaming" mode, -2L for "strict subsetting" mode (no new
### seqlevels added), -1L for "extended subsetting" mode (new seqlevels added),
### or an integer vector containing the mapping from the new to the old
### seqlevels for "general" mode (i.e. a combination of renaming and/or
### subsetting). Note that the vector describing the "general" mode is
### guaranteed to contain no negative values.
getSeqlevelsReplacementMode <- function(new_seqlevels, old_seqlevels)
{
    if (!is.character(new_seqlevels)
     || anyNA(new_seqlevels)
     || anyDuplicated(new_seqlevels))
        stop(wmsg("the supplied 'seqlevels' must be a character vector ",
                  "with no NAs and no duplicates"))
    nsl_names <- names(new_seqlevels)
    if (!is.null(nsl_names)) {
        nonempty_names <- nsl_names[!(nsl_names %in% c(NA, ""))]
        if (any(duplicated(nonempty_names)) ||
            length(setdiff(nonempty_names, old_seqlevels)) != 0L)
            stop(wmsg("the names of the supplied 'seqlevels' contain ",
                      "duplicates or invalid sequence levels"))
        return(match(nsl_names, old_seqlevels))
    }
    if (all(new_seqlevels %in% old_seqlevels))
        return(-2L)
    if (length(new_seqlevels) != length(old_seqlevels))
        return(-1L)
    is_renamed <- new_seqlevels != old_seqlevels
    tmp <- intersect(new_seqlevels[is_renamed], old_seqlevels[is_renamed])
    if (length(tmp) != 0L)
        return(-1L)
    return(-3L)
}

### Returns a logical vector of the same length as 'new_seqinfo' indicating
### whether the length or circularity flag of the corresponding sequence
### have changed. Assumes that the 'new2old' mapping is valid (see
### .reverse_new2old() function above for what this means exactly).
### NAs in 'new2old' are propagated to the result.
sequenceGeometryHasChanged <- function(new_seqinfo, old_seqinfo, new2old=NULL)
{
    ans_len <- length(new_seqinfo)
    if (is.null(new2old)) {
        idx1 <- idx2 <- seq_len(ans_len)
    } else {
        idx1 <- which(!is.na(new2old))
        idx2 <- new2old[idx1]
    }
    new_seqlengths <- seqlengths(new_seqinfo)[idx1]
    old_seqlengths <- seqlengths(old_seqinfo)[idx2]
    new_isCircular <- isCircular(new_seqinfo)[idx1]
    old_isCircular <- isCircular(old_seqinfo)[idx2]
    hasNotChanged <- function(x, y)
        (is.na(x) & is.na(y)) | (is.na(x) == is.na(y) & x == y)
    ans <- logical(ans_len)
    ans[] <- NA
    ans[idx1] <- !(hasNotChanged(new_seqlengths, old_seqlengths) &
                   hasNotChanged(new_isCircular, old_isCircular))
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqinfo() getter and setter
###

setGeneric("seqinfo", function(x) standardGeneric("seqinfo"))

setGeneric("seqinfo<-", signature="x",
    function(x, new2old=NULL,
             pruning.mode=c("error", "coarse", "fine", "tidy"),
             value)
        standardGeneric("seqinfo<-")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqnames() getter and setter
###

setGeneric("seqnames", function(x) standardGeneric("seqnames"))

setGeneric("seqnames<-", signature="x",
    function(x, value) standardGeneric("seqnames<-")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlevels() getter and setter
###

setGeneric("seqlevels", function(x) standardGeneric("seqlevels"))

### Default "seqlevels" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("seqlevels", "ANY", function(x) seqlevels(seqinfo(x)))

setGeneric("seqlevels<-", signature="x",
    function(x,
             pruning.mode=c("error", "coarse", "fine", "tidy"),
             value)
        standardGeneric("seqlevels<-")
)

### Default "seqlevels<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("seqlevels", "ANY",
    function(x,
             pruning.mode=c("error", "coarse", "fine", "tidy"),
             value)
    {
        ## Make the new Seqinfo object.
        x_seqinfo <- seqinfo(x)
        seqlevels(x_seqinfo) <- value
        ## Map the new sequence levels to the old ones.
        new2old <- getSeqlevelsReplacementMode(value, seqlevels(x))
        if (identical(new2old, -3L)) {
            ## "renaming" mode
            new2old <- seq_along(value)
        } else if (identical(new2old, -2L) || identical(new2old, -1L)) {
            ## "subsetting" mode
            new2old <- match(value, seqlevels(x))
        }
        ## Do the replacement.
        seqinfo(x, new2old=new2old, pruning.mode=pruning.mode) <- x_seqinfo
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sortSeqlevels()
###

setGeneric("sortSeqlevels", signature="x",
    function(x, X.is.sexchrom=NA) standardGeneric("sortSeqlevels")
)

setMethod("sortSeqlevels", "character",
    function(x, X.is.sexchrom=NA)
    {
        x[order(rankSeqlevels(x, X.is.sexchrom=X.is.sexchrom))]
    }
)

setMethod("sortSeqlevels", "ANY",
    function(x, X.is.sexchrom=NA)
    {
        seqlevels(x) <- sortSeqlevels(seqlevels(x),
                                      X.is.sexchrom=X.is.sexchrom)
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlevelsInUse() getter
###

setGeneric("seqlevelsInUse", function(x) standardGeneric("seqlevelsInUse"))

### Covers GenomicRanges, SummarizedExperiment, GAlignments, and any object
### for which the seqnames are returned as a factor-Rle.
setMethod("seqlevelsInUse", "Vector",
    function(x)
    {
        f <- runValue(seqnames(x))
        levels(f)[tabulate(f, nbins=nlevels(f)) != 0L]
    }
)

### Covers GRangesList and GAlignmentsList objects.
setMethod("seqlevelsInUse", "CompressedList",
    function(x) seqlevelsInUse(x@unlistData)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlevels0() getter
###
### Currently applicable to TxDb objects only.
###

setGeneric("seqlevels0", function(x) standardGeneric("seqlevels0"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlengths() getter and setter
###

setGeneric("seqlengths", function(x) standardGeneric("seqlengths"))

### Default "seqlengths" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("seqlengths", "ANY", function(x) seqlengths(seqinfo(x)))

setGeneric("seqlengths<-", signature="x",
    function(x, value) standardGeneric("seqlengths<-")
)

### Default "seqlengths<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("seqlengths", "ANY",
    function(x, value)
    {
        seqlengths(seqinfo(x)) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isCircular() getter and setter
###

setGeneric("isCircular", function(x) standardGeneric("isCircular"))

### Default "isCircular" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("isCircular", "ANY", function(x) isCircular(seqinfo(x)))

setGeneric("isCircular<-", signature="x",
    function(x, value) standardGeneric("isCircular<-")
)

### Default "isCircular<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("isCircular", "ANY",
    function(x, value)
    {
        isCircular(seqinfo(x)) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### genome() getter and setter
###

setGeneric("genome", function(x) standardGeneric("genome"))

### Default "genome" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("genome", "ANY", function(x) genome(seqinfo(x)))

setGeneric("genome<-", signature="x",
    function(x, value) standardGeneric("genome<-")
)

### Default "genome<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("genome", "ANY",
    function(x, value)
    {
        genome(seqinfo(x)) <- value
        x
    }
)

