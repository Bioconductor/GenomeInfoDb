### =========================================================================
### The "GenomeDescription" class
### -------------------------------------------------------------------------

setClass("GenomeDescription",
    representation(
        ## organism: "Homo sapiens", "Mus musculus", etc...
        organism="character",

        ## species: "Human", "Mouse", etc...
        species="character",

        ## provider: "UCSC", "BDGP", etc...
        provider="character",

        ## provider_version: "hg18", "mm8", "sacCer1", etc...
        provider_version="character",

        ## release_date: "Mar. 2006", "Feb. 2006", "Oct. 2003", etc...
        release_date="character",

        ## release_name: "NCBI Build 36.1", "NCBI Build 36",
        ## "SGD 1 Oct 2003 sequence", etc...
        release_name="character",

        ## names, lengths, and circularity flags of the genome sequences
        seqinfo="Seqinfo"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("organism", function(x) standardGeneric("organism"))
setMethod("organism", "GenomeDescription", function(x) x@organism)

setGeneric("species", function(x) standardGeneric("species"))
setMethod("species", "GenomeDescription", function(x) x@species)

setGeneric("provider", function(x) standardGeneric("provider"))
setMethod("provider", "GenomeDescription", function(x) x@provider)

setGeneric("providerVersion", function(x) standardGeneric("providerVersion"))
setMethod("providerVersion", "GenomeDescription", function(x) x@provider_version)

setGeneric("releaseDate", function(x) standardGeneric("releaseDate"))
setMethod("releaseDate", "GenomeDescription", function(x) x@release_date)

setGeneric("releaseName", function(x) standardGeneric("releaseName"))
setMethod("releaseName", "GenomeDescription", function(x) x@release_name)

setGeneric("bsgenomeName", function(x) standardGeneric("bsgenomeName"))
setMethod("bsgenomeName", "GenomeDescription",
    function(x)
    {
        part1 <- "BSgenome"
        tmp <- strsplit(organism(x), " ", fixed=TRUE)[[1L]]
        part2 <- paste(substr(tmp[1L], start=1L, stop=1L), tmp[2L], sep="")
        part3 <- provider(x)
        part4 <- providerVersion(x)
        paste(part1, part2, part3, part4, sep=".")
    }
)

setMethod("seqinfo", "GenomeDescription", function(x) x@seqinfo)

setMethod("seqnames", "GenomeDescription",
    function(x)
    {
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

GenomeDescription <- function(organism, species,
                              provider, provider_version,
                              release_date, release_name,
                              seqinfo)
{
    if (identical(organism, "NA")) organism <- NA_character_
    if (identical(species, "NA")) species <- NA_character_
    if (identical(release_date, "NA")) release_date <- NA_character_
    if (identical(release_name, "NA")) release_name <- NA_character_
    new("GenomeDescription",
        organism=organism,
        species=species,
        provider=provider,
        provider_version=provider_version,
        release_date=release_date,
        release_name=release_name,
        seqinfo=seqinfo)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method
###

### Kind of very low-level. Could go into S4Vectors if someone else needed
### this...
.compactPrintNamedAtomicVector <- function(x, margin="")
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
    cat(BiocGenerics:::labeledLine(margin, showMatrix[1L, ], count=FALSE,
                                           labelSep=""), sep="")
    cat(BiocGenerics:::labeledLine(margin, showMatrix[2L, ], count=FALSE,
                                           labelSep=""), sep="")
}

### NOT exported (but used in the BSgenome package).
showGenomeDescription <- function(x, margin="", print.seqlengths=FALSE)
{
    cat(margin, "organism: ", organism(x), " (",  species(x), ")\n", sep="")
    cat(margin, "provider: ", provider(x), "\n", sep="")
    cat(margin, "provider version: ", providerVersion(x), "\n", sep="")
    cat(margin, "release date: ", releaseDate(x), "\n", sep="")
    cat(margin, "release name: ", releaseName(x), "\n", sep="")
    if (print.seqlengths) {
        cat(margin, "---\n", sep="")
        cat(margin, "seqlengths:\n", sep="")
        .compactPrintNamedAtomicVector(seqlengths(x), margin=margin)
    }
}

setMethod("show", "GenomeDescription",
    function(object)
    {
        showGenomeDescription(object, margin="| ", print.seqlengths=TRUE)
    }
)

