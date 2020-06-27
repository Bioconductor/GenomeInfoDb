### =========================================================================
### The "GenomeDescription" class
### -------------------------------------------------------------------------

setClass("GenomeDescription",
    representation(
        ## organism: "Homo sapiens", "Mus musculus", etc...
        organism="character",

        ## common_name: "Human", "Mouse", etc...
        common_name="character",

        ## provider: "UCSC", "BDGP", etc...
        provider="character",

        ## release_date: "Mar. 2006", "Feb. 2006", "Oct. 2003", etc...
        release_date="character",

        ## names, lengths, and circularity flags of the genome sequences
        seqinfo="Seqinfo"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("organism", "GenomeDescription", function(object) object@organism)

setGeneric("commonName", function(object) standardGeneric("commonName"))
setMethod("commonName", "GenomeDescription",
    function(object) object@common_name
)

setGeneric("provider", function(x) standardGeneric("provider"))
setMethod("provider", "GenomeDescription", function(x) x@provider)

setGeneric("providerVersion", function(x) standardGeneric("providerVersion"))
setMethod("providerVersion", "GenomeDescription",
    function(x)
    {
        .Deprecated()
        genome(x)[[1L]]
    }
)

setGeneric("releaseDate", function(x) standardGeneric("releaseDate"))
setMethod("releaseDate", "GenomeDescription", function(x) x@release_date)

setGeneric("bsgenomeName", function(x) standardGeneric("bsgenomeName"))
setMethod("bsgenomeName", "GenomeDescription",
    function(x)
    {
        part1 <- "BSgenome"
        tmp <- strsplit(organism(x), " ", fixed=TRUE)[[1L]]
        part2 <- paste(substr(tmp[1L], start=1L, stop=1L), tmp[2L], sep="")
        part3 <- provider(x)
        part4 <- genome(x)[[1L]]
        paste(part1, part2, part3, part4, sep=".")
    }
)

setMethod("seqinfo", "GenomeDescription", function(x) x@seqinfo)

setMethod("seqnames", "GenomeDescription",
    function(x)
    {
        ## Do NOT use 'seqnames(x)' here or you'll get infinite recursion!
        ans <- seqnames(seqinfo(x))
        if (length(ans) == 0L)
            ans <- NULL
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

setValidity("GenomeDescription",
    function(object)
    {
        SINGLE_STRING_SLOTS <- setdiff(slotNames("GenomeDescription"),
                                       "seqinfo")
        .validSlot <- function(slotname)
        {
            slotval <- slot(object, slotname)
            if (isSingleStringOrNA(slotval))
                return(NULL)
            problem <- paste("slot '", slotname, "' must be a ",
                             "single string (or NA)", sep="")
            return(problem)
        }
        unlist(lapply(SINGLE_STRING_SLOTS, .validSlot))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor-like functions
###

### NOTE: In BioC 3.1, the 'species' argument was replaced with the
### 'common_name' argument but the former was kept for backward compatibility
### (essentially with existing SNPlocs and XtraSNPlocs packages).
### TODO: At some point the 'species' argument needs to be deprecated.
GenomeDescription <- function(organism, common_name, provider, release_date,
                              seqinfo, species=NA_character_)
{
    if (identical(organism, "NA")) organism <- NA_character_
    if (missing(common_name))
        common_name <- species
    if (identical(common_name, "NA")) common_name <- NA_character_
    if (identical(release_date, "NA")) release_date <- NA_character_

    ## Check 'seqinfo'.
    if (!is(seqinfo, "Seqinfo"))
        stop(wmsg("'seqinfo' must be a Seqinfo object"))
    ugenome <- unique(genome(seqinfo))
    if (length(ugenome) != 1L || is.na(ugenome))
        stop(wmsg("'genome(seqinfo)' must contain a unique string"))

    new("GenomeDescription",
        organism=organism,
        common_name=common_name,
        provider=provider,
        release_date=release_date,
        seqinfo=seqinfo)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method
###

### NOT exported but used in BSgenome package.
### Kind of very low-level. Could go into S4Vectors if someone else needed
### this...
compactPrintNamedAtomicVector <- function(x, margin="")
{
    x_len <- length(x)
    halfWidth <- (getOption("width") - nchar(margin)) %/% 2L
    first <- max(1L, halfWidth)
    showMatrix <-
      rbind(as.character(head(names(x), first)),
            as.character(head(x, first)))
    if (x_len > first) {
        last <- min(x_len - first, halfWidth)
        showMatrix <-
          cbind(showMatrix,
                rbind(as.character(tail(names(x), last)),
                      as.character(tail(x, last))))
    }
    showMatrix <- format(showMatrix, justify="right")
    cat(S4Vectors:::labeledLine(margin, showMatrix[1L, ], count=FALSE,
                                        labelSep=""), sep="")
    cat(S4Vectors:::labeledLine(margin, showMatrix[2L, ], count=FALSE,
                                        labelSep=""), sep="")
}

### NOT exported (but used in the BSgenome package).
showGenomeDescription <- function(x, margin="", print.seqlengths=FALSE)
{
    cat(margin, "organism: ", organism(x), " (",  commonName(x), ")\n", sep="")
    cat(margin, "genome: ", genome(x)[[1L]], "\n", sep="")
    cat(margin, "provider: ", provider(x), "\n", sep="")
    cat(margin, "release date: ", releaseDate(x), "\n", sep="")
    if (print.seqlengths) {
        cat(margin, "---\n", sep="")
        cat(margin, "seqlengths:\n", sep="")
        compactPrintNamedAtomicVector(seqlengths(x), margin=margin)
    }
}

setMethod("show", "GenomeDescription",
    function(object)
    {
        showGenomeDescription(object, margin="| ", print.seqlengths=TRUE)
    }
)

