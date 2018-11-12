test_listOrganisms <- function(){

    # check the file exists
    filename <- system.file(package="GenomeInfoDb",  "extdata",
                            "dataFiles", "genomeMappingTbl.csv")
    checkTrue(file.exists(filename))

    # check format of file
    data <- read.csv(filename, header=TRUE, stringsAsFactors=FALSE)
    checkIdentical(c("commonName", "organism", "ensemblID",
                     "origEnsemblVersion", "origEnsemblDate", "ensemblVersion",
                     "ensemblDate", "ucscID", "ucscDate", "releaseName"),
                   colnames(data))

    # check helper function listOrganisms
    checkIdentical(class(listOrganisms()), "data.frame")
    checkTrue(all(dim(listOrganisms()) != c(0,0)))

}

filename <- system.file(package="GenomeInfoDb",  "extdata",
                        "dataFiles", "genomeMappingTbl.csv")
data <- read.csv(filename, header=TRUE, stringsAsFactors=FALSE)

test_genomeBuilds <- function(){

    # returns ERROR species not identified
    checkException(genomeBuilds(style="UCSC"))

    # returns default style UCSC for mouse
    testBuild <- genomeBuilds("Mouse")
    checkIdentical(unique(testBuild$commonName), "mouse")
    tblid <- sort(unique(tolower(data$ucscID)[which(
        tolower(data$organism)==tolower("Mus musculus"))]))
    checkIdentical(tblid, sort(tolower(unique(testBuild$ucscID))))

    # test multiple spcies and specify ucsc
    testBuild <- genomeBuilds(c("Mouse", "Dog"), style="UCSC")
    checkIdentical(unique(testBuild$commonName),c("dog","mouse"))
    tblid <- sort(unique(tolower(data$ucscID)[which(
        tolower(data$organism) ==tolower("Mus musculus") |
        tolower(data$commonName) =="dog")]))
    checkIdentical(tblid, sort(tolower(unique(testBuild$ucscID))))

    # test Ensembl style
    testBuild <- genomeBuilds(c("Mouse", "Dog"), style="Ensembl")
    checkIdentical(unique(testBuild$commonName),c("dog","mouse"))
    tblid <- sort(unique(tolower(data$ensemblID)[which(
        tolower(data$commonName) =="mouse" |
        tolower(data$commonName) =="dog")]))
    checkIdentical(tblid, sort(tolower(unique(testBuild$ensemblID))))

    # test species not found
    testBuild <- suppressWarnings({
        genomeBuilds(c("Mouse", "NotHere"), style="Ensembl")
    })
    checkIdentical(unique(testBuild$commonName), "mouse")
    warn <- tryCatch(genomeBuilds(c("Mouse", "NotHere"), style="Ensembl"),
                     warning=conditionMessage)
    checkIdentical(warn, "'organism' not found: nothere")
    testBuild <- suppressWarnings(genomeBuilds("NotHere"))
    checkTrue(all(dim(testBuild) == c(0L, 0L)))
    warn <- tryCatch(genomeBuilds("NotHere"), warning=conditionMessage)
    checkIdentical(warn, "'organism' not found: nothere")

    # test style not found - ERROR
    checkException(genomeBuilds("Mouse", style="NotHere"), silent=TRUE)

}



test_mapGenomeBuilds <- function(){

    # returns ERROR genome not identified
    checkException(mapGenomeBuilds())

    # returns default style UCSC
    testBuild <- mapGenomeBuilds("NCBIm37")
    checkIdentical(unique(testBuild$commonName), "mouse")
    idx = c(which(data$ensemblID == "NCBIm37"), which(data$ucscID == "NCBIm37"))
    checkIdentical(sort(unique(data$ucscID[idx])),
                   sort(unique(testBuild$ucscID)))

    # test multiple genome and specify ucsc
    testBuild <- mapGenomeBuilds(c("NCBIm37", "canFam3"), style="UCSC")
    checkIdentical(unique(testBuild$commonName), c("dog", "mouse"))
    buildVal <- sort(unique(testBuild$ucscID))
    idx = c(which(data$ensemblID == "NCBIm37"), which(data$ucscID == "NCBIm37"),
        which(data$ensemblID == "canFam3"), which(data$ucscID == "canFam3"))
    checkIdentical(sort(unique(data$ucscID[idx])), buildVal)

    # test Ensembl style
    testBuild <- mapGenomeBuilds(c("NCBIm37", "canFam3"), style="Ensembl")
    checkIdentical(unique(testBuild$organism),
                   c("Canis familiaris", "Mus musculus"))
    buildVal <- sort(unique(unique(testBuild$ensemblID)))
    idx = c(which(data$ensemblID == "NCBIm37"), which(data$ucscID == "NCBIm37"),
        which(data$ensemblID == "canFam3"), which(data$ucscID == "canFam3"))
    checkIdentical(sort(unique(data$ensemblID[idx])), buildVal)

    # test genome not found
    testBuild <- suppressWarnings({
        mapGenomeBuilds(c("canFam3", "NotHere"), style="Ensembl")
    })
        
    checkIdentical(unique(testBuild$commonName), "dog")
    warn <- tryCatch(mapGenomeBuilds(c("canFam3", "NotHere"), style="Ensembl"),
                     warning=conditionMessage)
    checkIdentical(warn, "'genome' not found: nothere")
    testBuild <- suppressWarnings(mapGenomeBuilds("NotHere"))
    checkTrue(all(dim(testBuild) == c(0L, 0L)))
    warn <- tryCatch(mapGenomeBuilds("NotHere"), warning=conditionMessage)
    checkIdentical(warn, "'genome' not found: nothere")

    # test style not found - ERROR
    checkException(mapGenomeBuilds("NCBIm37", style="NotHere"))
}
