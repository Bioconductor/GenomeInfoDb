.getDatadir <-
    function()
{
    system.file(package = "GenomeInfoDb","extdata","dataFiles")    
}        

.getNamedFiles <-
    function()
{
    filePath <- .getDatadir()
    files <- dir(filePath, full.names=TRUE)
    setNames(files, sub(".txt$", "", basename(files)))
}

.getDataInFile <- 
    function(name)
{
    ##name will always be informat Homo sapiens
    ## files are in format of : Homo_sapiens.txt
    filename <- paste0(.getDatadir(),"/",sub(" ", "_", name),".txt")
    if (file.exists(filename)) {
        read.table(filename, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    } else {
        stop("The species, " ,name, " is not supported by GenomeInfoDb")
    }
    
}  

.isSingleString <- function (x)
{
    is.character(x) && length(x) == 1L && !is.na(x)
}

.isTRUEorFALSE <- function (x) 
{
    is.logical(x) && length(x) == 1L && !is.na(x)
}


.supportedSeqnameStyles <- 
    function()
{
    dom <- lapply(.getNamedFiles(), scan, nlines=1, what=character(),
                  quiet=TRUE)
    lapply(dom, function(x) {x[!(x %in% c("linear","auto","sex"))] })
}


.isSupportedSeqnamesStyle <- 
    function(species, style)
{
    if (missing(species) || missing(style))
        stop("'species' or 'style' missing")
    
    species <- sub(" ", "_", species)
    possible <- lapply(.getNamedFiles(), scan, nlines=1, what=character(),
                       quiet=TRUE)
    availStyles <- possible[[species]]
    style %in% availStyles
}

.isSupportedSeqnames <- 
    function(species, style, seqnames)
{
    if (missing(species) || missing(style) || missing(seqnames))
        stop("'species', 'style' or 'seqnames' missing")    
    
    trueSeq <- extractSeqlevels(species=species,style=style)
    all(trueSeq %in% seqnames)
}

.supportedSeqnameMappings <-
    function()
{
    dom <-  lapply(.getNamedFiles(), read.table, header=TRUE, sep="\t",
                   stringsAsFactors=FALSE)
    lapply(dom, function(x) {x[,-c(1:3)] })
}

.guessSpeciesStyle <- 
    function(seqnames)
{
    zz <- .supportedSeqnameMappings()
    got2 <- lapply(zz ,function(y) lapply(y, function(z) 
        sum(z %in% seqnames)) )
    unlistgot2 <- unlist(got2, recursive=TRUE,use.names=TRUE)
    
    if (max(unlistgot2) == 0) {
        txt <- "The style does not have a compatible entry for the
    species supported by Seqname. Please see
    genomeStyles() for supported species/style"
        stop(paste(strwrap(txt, exdent=2), collapse="\n"))
    }
    
    ##vec is in format "Homo_sapiens.UCSC"
    vec <- names(which.max(unlistgot2))  
    species <- sub("_", " ",unlist(strsplit(vec,"[.]")),fixed=TRUE)[[1]]
    style <- unlist(strsplit(vec,"[.]"))[[2]]      
    c(species,style)
}        


