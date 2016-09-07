listSpecies <- function(){

    filename <- system.file(package="GenomeInfoDb",  "extdata",
                            "dataFiles", "genomeMappingTbl.csv")
    tbl <- read.csv(filename, header=TRUE, stringsAsFactors=FALSE)
    sort(unique(tbl$speciesShort))
}

genomeBuilds <- function(species, style = c("UCSC", "Ensembl")) {

    if (!is.character(species))
        stop("'species' must be a character vector")
    style <- match.arg(style)

    filename <- system.file(package="GenomeInfoDb",  "extdata",
                            "dataFiles", "genomeMappingTbl.csv")
    tbl <- read.csv(filename, header=TRUE, stringsAsFactors=FALSE)

    colkeep <- switch(style,
                   UCSC="ucscID",
                   Ensembl="ensemblID"
               )

    fnd <- sapply(tolower(species), grep, tolower(tbl$speciesShort))
    notFnd <- names(which(lengths(fnd) == 0))
    if (length(notFnd))
        warning("'species' not found: ", paste(notFnd, collapse=", "),
                call.=FALSE)

    if (!missing(species))
        tbl <- tbl[tolower(tbl$speciesShort) %in% tolower(species),
                   c("speciesShort", colkeep)]
    else
        tbl <- tbl[, c("speciesShort", colkeep)]
    if (nrow(tbl) == 0L)
        return(list())

    rownames(tbl) <- NULL
    lst <- lapply(unique(tolower(tbl$speciesShort)),
                  function(xx) unique(na.omit(tbl[tolower(tbl$speciesShort) ==
                                                  xx, colkeep])))
    names(lst) <- unique(tolower(tbl$speciesShort))
    lst
}

mapGenomeBuilds <- function(genome, style = c("UCSC", "Ensembl") ){

    if (!is.character(genome))
            stop("'genome' must be a character vector")
    genome <- tolower(genome)
    style <- match.arg(style)

    filename <- system.file(package="GenomeInfoDb",  "extdata",
                            "dataFiles", "genomeMappingTbl.csv")
    tbl <- read.csv(filename, header=TRUE, stringsAsFactors=FALSE)

    colkeep <- switch(tolower(style),
                      ucsc=c("ucscID","ucscDate","ensemblID"),
                      ensembl=c("ensemblID","ensemblVersion", "ensemblDate",
                          "ucscID" ))

    fnd = sapply(genome, grep, tolower(c(tbl$ensemblID, tbl$ucscID)))
    notFnd = names(which(sapply(FUN=length, fnd)==0))
    if(length(notFnd) != 0L)
        warning("'genome' not found: ", paste(notFnd,collapse=", "),
                call.=FALSE)

    rowkeep <- ((tolower(tbl$ucscID) %in% genome) |
                (tolower(tbl$ensemblID) %in% genome))
    if (sum(rowkeep) == 0)
        return(list())
    tbl <- tbl[rowkeep, c("speciesShort", colkeep)]

    rownames(tbl) <- NULL
    species = unique(tbl$speciesShort)
    tbl = lapply(species,
        function(xx) unique(tbl[tbl$speciesShort == xx,colkeep]))
    names(tbl) = species
    tbl
}
