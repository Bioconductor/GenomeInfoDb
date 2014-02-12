test_basic <- 
    function()
{
    filename <- system.file(package="Seqnames",  "extdata",
                            "dataFiles", "Homo_sapiens.txt")
    checkTrue(file.exists(filename))
   
    # check the format of the file
    data<- read.table(filename,header=TRUE,sep="\t")
    checkEquals(c(25,5),dim(data))
    checkEquals(c('linear','auto','sex','NCBI','UCSC'),names(data))
    
    #check if first 3 columns contain only true or false entries
    checkEquals(c(TRUE,FALSE),unique(data[,1]))
    checkEquals(c(TRUE,FALSE),unique(data[,2]))
    checkEquals(c(FALSE,TRUE),unique(data[,3]))
}

test_supportedStyles <- 
    function()
{
    checkEquals("data.frame", class(supportedStyles("Homo sapiens")))
    checkEquals(c(25,5), dim(supportedStyles("Homo sapiens")))
    checkException(supportedStyles("SAD"))
}

test_seqnamesOrder <- 
    function()
{
    checkIdentical(integer(), seqnamesOrder())
    checkException(seqnamesOrder(c(1,2,3,4,5)))
    
    seqnames <- c("chr1","chr9", "chr2", "chr3", "chr10")
    checkEquals(c(1L, 3L, 4L, 2L, 5L), seqnamesOrder(seqnames) )
}    

test_extractSeqnameSet <- 
    function()
{
    got <- extractSeqnameSet("Homo sapiens", "UCSC" )
    checkEquals(25,length(got))
    checkEquals("character",class(got))
    
    checkException(extractSeqnameSet("aaa","Homo sapiens"))
    checkException(extractSeqnameSet("Drosophila melanogaster"))
}

test_extractSeqnameSetByGroup <- 
    function()
{
    got <- extractSeqnameSetByGroup("Drosophila melanogaster","Ensembl","auto")
    checkEquals(5,length(got))
    checkEquals("character",class(got))
    
    checkException(extractSeqnameSetByGroup("aaa","Homo sapiens"))
    checkException(extractSeqnameSetByGroup("Drosophila melanogaster"))
}

test_seqnamesInGroup <-
    function()
{
    newch <- paste0("chr",c(1:22,"X","Y","M","1_gl000192_random","4_ctg9_hap1"))
    got1 <- seqnamesInGroup(newch, group="sex")
    checkEquals(c("chrX","chrY"),got1)
    
    newchr <- as.character(c(1:22,"X","Y","MT"))
    got2 <- seqnamesInGroup(newchr, group="all","Homo sapiens","NCBI")
    checkEquals(25,length(got2))
}  

test_seqnameStyle <-
    function()
{
       checkEquals("NCBI",seqnameStyle(c("2L","2R","X","Xhet")))
       checkEquals("UCSC",seqnameStyle(paste0("chr",c(1:30))))
}
