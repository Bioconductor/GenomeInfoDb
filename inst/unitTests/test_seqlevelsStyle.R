test_basic <- 
    function()
{
    filename <- system.file(package="GenomeInfoDb",  "extdata",
                            "dataFiles", "Homo_sapiens.txt")
    checkTrue(file.exists(filename))
   
    # check the format of the file
    data <- read.table(filename,header=TRUE,sep="\t")
    checkIdentical(c(25L, 6L), dim(data))
    checkIdentical(c('circular', 'auto', 'sex', 'NCBI', 'UCSC', 'dbSNP'),
                   colnames(data))
    
    #check if first 3 columns contain only true or false entries
    checkEquals(c(FALSE,TRUE),unique(data[,1]))
    checkEquals(c(TRUE,FALSE),unique(data[,2]))
    checkEquals(c(FALSE,TRUE),unique(data[,3]))
}

test_genomeStyles <- 
    function()
{
    checkIdentical("data.frame", class(genomeStyles("Homo sapiens")))
    checkIdentical(c(25L, 6L), dim(genomeStyles("Homo sapiens")))
    checkException(genomeStyles("SAD"))
}

test_orderSeqlevels <- 
    function()
{
    checkIdentical(integer(), orderSeqlevels())
    checkException(orderSeqlevels(c(1,2,3,4,5)))
    
    seqnames <- c("chr1","chr9", "chr2", "chr3", "chr10")
    checkEquals(c(1L, 3L, 4L, 2L, 5L), orderSeqlevels(seqnames) )
}    

test_extractSeqlevels <- 
    function()
{
    got <- extractSeqlevels("Homo sapiens", "UCSC" )
    checkEquals(25,length(got))
    checkEquals("character",class(got))
    
    checkException(extractSeqlevels("aaa","Homo sapiens"))
    checkException(extractSeqlevels("Drosophila melanogaster"))
}

test_extractSeqlevelsByGroup <- 
    function()
{
    got <- extractSeqlevelsByGroup("Drosophila melanogaster","Ensembl","auto")
    checkEquals(5,length(got))
    checkEquals("character",class(got))
    
    checkException(extractSeqlevelsByGroup("aaa","Homo sapiens"))
    checkException(extractSeqlevelsByGroup("Drosophila melanogaster"))
    checkException(extractSeqlevelsByGroup("Homo sapiens","auto","NCBI"))
}

test_seqlevelsInGroup <-
    function()
{
    newch <- paste0("chr",c(1:22,"X","Y","M","1_gl000192_random","4_ctg9_hap1"))
    got1 <- seqlevelsInGroup(newch, group="sex")
    checkEquals(c("chrX","chrY"),got1)
    
    newchr <- as.character(c(1:22,"X","Y","MT"))
    got2 <- seqlevelsInGroup(newchr, group="all","Homo sapiens","NCBI")
    checkEquals(25,length(got2))
}  

test_seqlevelsStyle <-
    function()
{
       checkEquals("NCBI",seqlevelsStyle(c("2L","2R","X","Xhet")))
       checkEquals("UCSC",seqlevelsStyle(paste0("chr",c(1:30))))
}
