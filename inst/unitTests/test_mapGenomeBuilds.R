test_listSpecies <- function(){

    # check the file exists
    filename <- system.file(package="GenomeInfoDb",  "extdata",
                            "dataFiles", "genomeMappingTbl.csv")
    checkTrue(file.exists(filename))

    # check format of file
    data <- read.csv(filename, header=TRUE, stringsAsFactors=FALSE)
    checkIdentical(c("speciesShort", "speciesLong", "ensemblID",
                     "origEnsemblVersion", "origEnsemblDate", "ensemblVersion",
                     "ensemblDate", "ucscID", "ucscDate", "releaseName"),
                   colnames(data))

    # check helper function listSpecies
    checkIdentical(class(listSpecies()), "character")
    checkTrue(length(listSpecies()) > 0)

}

filename <- system.file(package="GenomeInfoDb",  "extdata",
                        "dataFiles", "genomeMappingTbl.csv")
data <- read.csv(filename, header=TRUE, stringsAsFactors=FALSE)

test_genomeBuilds <- function(){

    # returns ERROR species not identified
    checkException(genomeBuilds(style="UCSC"))

    # returns default style UCSC for mouse
    testBuild <- genomeBuilds("Mouse")
    checkIdentical(names(testBuild), "mouse")
    tblid <- sort(unique(tolower(data$ucscID)[which(
        tolower(data$speciesShort)=="mouse")]))
    checkIdentical(tblid, sort(tolower(unlist(testBuild, use.names=FALSE))))

    # test multiple spcies and specify ucsc
    testBuild <- genomeBuilds(c("Mouse", "Dog"), style="UCSC")
    checkIdentical(names(testBuild),c("dog","mouse"))
    tblid <- sort(unique(tolower(data$ucscID)[which(
        tolower(data$speciesShort) =="mouse" |
        tolower(data$speciesShort) =="dog")]))
    checkIdentical(tblid, sort(tolower(unlist(testBuild, use.names=FALSE))))

    # test Ensembl style
    testBuild <- genomeBuilds(c("Mouse", "Dog"), style="Ensembl")
    checkIdentical(names(testBuild),c("dog","mouse"))
    tblid <- sort(unique(tolower(data$ensemblID)[which(
        tolower(data$speciesShort) =="mouse" |
        tolower(data$speciesShort) =="dog")]))
    checkIdentical(tblid, sort(tolower(unlist(testBuild, use.names=FALSE))))

    # test species not found
    testBuild <- genomeBuilds(c("Mouse", "NotHere"), style="Ensembl")
    checkIdentical(names(testBuild), "mouse")
    warn <- tryCatch(genomeBuilds(c("Mouse", "NotHere"), style="Ensembl"),
                     warning=conditionMessage)
    checkIdentical(warn, "'species' not found: nothere")
    testBuild <- genomeBuilds("NotHere")
    checkTrue(length(testBuild) == 0L)
    warn <- tryCatch(genomeBuilds("NotHere"), warning=conditionMessage)
    checkIdentical(warn, "'species' not found: nothere")

    # test style not found - ERROR
    checkException(genomeBuilds("Mouse", style="NotHere"))

}



test_mapGenomeBuilds <- function(){

    # returns ERROR genome not identified
    checkException(mapGenomeBuilds())

    # returns default style UCSC
    testBuild <- mapGenomeBuilds("NCBIm37")
    checkIdentical(names(testBuild), "Mouse")
    idx = c(which(data$ensemblID == "NCBIm37"), which(data$ucscID == "NCBIm37"))
    checkIdentical(sort(unique(data$ucscID[idx])),
                   sort(unique(testBuild[[1]]$ucscID)))

    # test multiple genome and specify ucsc
    testBuild <- mapGenomeBuilds(c("NCBIm37", "canFam3"), style="UCSC")
    checkIdentical(names(testBuild), c("Dog", "Mouse"))
    buildVal <- sort(unique(unlist(testBuild, use.names=FALSE)[grepl(
                                                  names(unlist(testBuild)),
                                                  pattern="ucscID")]))
    idx = c(which(data$ensemblID == "NCBIm37"), which(data$ucscID == "NCBIm37"),
        which(data$ensemblID == "canFam3"), which(data$ucscID == "canFam3"))
    checkIdentical(sort(unique(data$ucscID[idx])), buildVal)

    # test Ensembl style
    testBuild <- mapGenomeBuilds(c("NCBIm37", "canFam3"), style="Ensembl")
    checkIdentical(names(testBuild), c("Dog", "Mouse"))
    buildVal <- sort(unique(unlist(testBuild, use.names=FALSE)[grepl(
                                                  names(unlist(testBuild)),
                                                  pattern="ensemblID")]))
    idx = c(which(data$ensemblID == "NCBIm37"), which(data$ucscID == "NCBIm37"),
        which(data$ensemblID == "canFam3"), which(data$ucscID == "canFam3"))
    checkIdentical(sort(unique(data$ensemblID[idx])), buildVal)

    # test genome not found
    testBuild <- mapGenomeBuilds(c("canFam3", "NotHere"), style="Ensembl")
    checkIdentical(names(testBuild), "Dog")
    warn <- tryCatch(mapGenomeBuilds(c("canFam3", "NotHere"), style="Ensembl"),
                     warning=conditionMessage)
    checkIdentical(warn, "'genome' not found: nothere")
    testBuild <- mapGenomeBuilds("NotHere")
    checkTrue(length(testBuild) == 0L)
    warn <- tryCatch(mapGenomeBuilds("NotHere"), warning=conditionMessage)
    checkIdentical(warn, "'genome' not found: nothere")

    # test style not found - ERROR
    checkException(mapGenomeBuilds("NCBIm37", style="NotHere"))

}
