supportedStyles <-
    function(species)
{
    if (missing(species))
        lapply(.getNamedFiles(), read.table, header=TRUE, sep="\t",
           stringsAsFactors=FALSE)
    else 
        .getDataInFile(species)
}


extractSeqnameSet <- 
    function(species, style)
{
    if (missing(species) || missing(style))
        stop("'species' or 'style' missing")    
    
    if(.isSupportedSeqnamesStyle(species, style))
    {
        data <- .getDataInFile(species)
        result <- as.vector(data[,which( names(data) %in% style)])
        
    }else{
        stop("The style specified by '",style,
             "' does not have a compatible entry for the species ",species)}   
    result
}

extractSeqnameSetByGroup <- 
    function(species, style, group)
{
    if (missing(species) || missing(style) || missing(group))
        stop("'species', 'style', and / or 'group' missing")     
        
    if(.isSupportedSeqnamesStyle(species, style))
    {
        data <- .getDataInFile(species)
        if (group!="all"){
            colInd <- which(names(data)%in% group)
            Ind <- which(data[,colInd]==1)
            result <- as.vector(data[Ind,which( names(data) %in% style)])
        }
        else{
            result <- as.vector(data[,which( names(data) %in% style)])
        }
    }else{
        stop("The style specified by '",style,
             "' does not have a compatible entry for the species ",species)}   
    result
}

findSequenceRenamingMaps <- 
    function(seqnames, style, best.only=TRUE, drop=TRUE)
{
    if (!is.character(seqnames))
        stop("'seqnames' must be a character vector")
    if (!.isSingleString(style))
        stop("the supplied seqname style must be a single string")
    if (!.isTRUEorFALSE(best.only))
        stop("'best.only' must be TRUE or FALSE")
    if (!.isTRUEorFALSE(drop))
        stop("'drop' must be TRUE or FALSE")
    supported_styles <- .supportedSeqnameStyles()
    tmp <- unlist(supported_styles, use.names = FALSE)
    compatible_species <- rep.int(names(supported_styles),
                                  sapply(supported_styles,NROW))
    compatible_species <- compatible_species[tolower(tmp) ==
                                                 tolower(style)]
    if (length(compatible_species) == 0L)
        stop("supplied seqname style \"", style, "\" is not supported")
    seqname_mappings <- .supportedSeqnameMappings()
    ans <- lapply(compatible_species, function(species) {
        mapping <- seqname_mappings[[species]]
        names(mapping) <- tolower(names(mapping))
        to_seqnames <- as.character(mapping[[tolower(style)]])
        lapply(mapping, function(from_seqnames) to_seqnames[match(seqnames,
                                                  from_seqnames)])
    })
    ans_ncol <- length(seqnames)
    ans <- matrix(unlist(ans, use.names = FALSE), ncol = ans_ncol,
                  byrow = TRUE)
    colnames(ans) <- seqnames
    score <- rowSums(!is.na(ans))
    idx <- score != 0L
    if (best.only)
        idx <- idx & (score == max(score))
    ans <- ans[idx, , drop = FALSE]
    ans <- as.matrix(unique(as.data.frame(ans, stringsAsFactors = FALSE)))
    if (nrow(ans) == 1L && drop)
        ans <- drop(ans)
    else rownames(ans) <- NULL
    ans        
}

seqnamesInGroup <- 
    function(seqnames, group=c("all", "auto", "sex", "linear"),
             species, style)
{
    group <- match.arg(group)
    if (missing(species) && missing(style)) {
        ## guess the species and / or style for the object
        ans <- .guessSpeciesStyle(seqnames)
        species<- ans[1]
        style <- ans[2]
    }
    
    if (.isSupportedSeqnamesStyle(species, style)) {
        seqvec <- extractSeqnameSetByGroup( species, style, group)
        seqvec[na.omit(match(seqnames, seqvec))]
    } else {
        txt <- paste0( "The style specified by ", sQuote(style),
                       " does not have a compatible entry for the species ",
                       sQuote(species))
        stop(paste(strwrap(txt, exdent=2), collapse="\n"))
    }        
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqnameStyle() getter and setter
###

setGeneric("seqnameStyle", 
    function(x) standardGeneric("seqnameStyle"))

### Default "seqnameStyle" method works on any object 'x' with a working
### "seqinfo" method; defined in GenomicRanges.

setMethod("seqnameStyle", "character",
    function(x) 
{
    ## implement seqnameStyle,character-method
    seqnames <- unique(x)      
    ans <- .guessSpeciesStyle(seqnames)
    ans[2]
})

setGeneric("seqnameStyle<-", signature="x",
    function(x, value) standardGeneric("seqnameStyle<-")
)

### Default "seqnameStyle<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods; defined in GenomicRanges.
