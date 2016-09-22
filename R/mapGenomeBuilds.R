listOrganisms <- function(){

    filename <- system.file(package="GenomeInfoDb",  "extdata",
                            "dataFiles", "genomeMappingTbl.csv")
    tbl <- read.csv(filename, header=TRUE, stringsAsFactors=FALSE)
    tbl_names <- unique(tbl[,1:2])
    rownames(tbl_names) <- NULL
    tbl_names[,2] = paste0(toupper(substring(tbl_names[,2], 1, 1)),
                 substring(tbl_names[,2], 2, nchar(tbl_names[,2])))
    tbl_names
}

genomeBuilds <- function(organism, style = c("UCSC", "Ensembl")) {

    if (!is.character(organism))
        stop("'organism' must be a character vector")
    style <- match.arg(style)

    filename <- system.file(package="GenomeInfoDb",  "extdata",
                            "dataFiles", "genomeMappingTbl.csv")
    tbl <- read.csv(filename, header=TRUE, stringsAsFactors=FALSE)

    colkeep <- switch(style,
                   UCSC="ucscID",
                   Ensembl="ensemblID"
               )

    fnd1 <- sapply(tolower(organism), grep, tolower(tbl$commonName))
    fnd2 <- sapply(tolower(organism), grep, tolower(tbl$organism))
    fnd <- mapply(c, fnd1, fnd2)
    notFnd <- names(which(lengths(fnd) == 0))
    if (length(notFnd))
        warning("'organism' not found: ", paste(notFnd, collapse=", "),
                call.=FALSE)

    if (!missing(organism)){
        rowkeep <- apply(FUN=any, MARGIN=1, cbind(
                         tolower(tbl$commonName) %in% tolower(organism),
                         tolower(tbl$organism) %in% tolower(organism)
                        ))             
        tbl <- tbl[rowkeep,c("commonName", "organism", colkeep)]
    }else
        tbl <- tbl[,c("commonName", "organism", colkeep)]

    if (nrow(tbl) == 0L)
        return(data.frame())

    tbl <- unique(na.omit(tbl))
    rownames(tbl) <- NULL
    tbl[,2] = paste0(toupper(substring(tbl[,2], 1, 1)),
                 substring(tbl[,2], 2, nchar(tbl[,2])))
    tbl
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
        return(data.frame())
    tbl <- tbl[rowkeep, c("commonName", "organism", colkeep)]

    rownames(tbl) <- NULL
    tbl[,2] = paste0(toupper(substring(tbl[,2], 1, 1)),
                 substring(tbl[,2], 2, nchar(tbl[,2])))
    tbl
}
