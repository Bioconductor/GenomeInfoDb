test_keepSeqlevels <- function()
{
    gr <- GRanges(c("chr1", "chr1", "chr2", "chr3"), IRanges(1:4, width=3))
    ## unnamed
    checkIdentical("chr1", seqlevels(keepSeqlevels(gr, "chr1")))
    got <- seqlevels(keepSeqlevels(gr, c("chr1", "chr3")))
    checkIdentical(c("chr1", "chr3"), got)
    ## named
    checkIdentical("chr1", seqlevels(keepSeqlevels(gr, c(foo="chr1"))))
    ## bogus
    got <- seqlevels(suppressWarnings(keepSeqlevels(gr, "chrX")))
    checkIdentical(character(0), got)
}

test_dropSeqlevels <- function()
{
    gr <- GRanges(c("chr1", "chr1", "chr2", "chr3"), IRanges(1:4, width=3))
    ## unnamed
    checkIdentical(c("chr2", "chr3"), seqlevels(dropSeqlevels(gr, "chr1")))
    got <- seqlevels(dropSeqlevels(gr, c("chr1", "chr3")))
    checkIdentical("chr2", got)
    ## named
    got <- seqlevels(dropSeqlevels(gr, c(foo="chr1")))
    checkIdentical(c("chr2", "chr3"), got)
    ## bogus
    got <- seqlevels(suppressWarnings(dropSeqlevels(gr, "chrX")))
    checkIdentical(seqlevels(gr), got)

    grl <- split(gr, as.character(seqnames(gr)))
    got <- dropSeqlevels(grl, c("chr1", "chr2"))
    checkIdentical("chr3", seqlevels(got))
    checkIdentical(1L, length(got))
}

test_renameSeqlevels <- function()
{
    opt <- options(useFancyQuotes=FALSE)
    on.exit(options(opt))

    gr <- GRanges(c("chr1", "chr1", "chr2", "chr3"), IRanges(1:4, width=3))
    checkException(renameSeqlevels(gr, "CHR1"), silent=TRUE)
    got <- seqlevels(renameSeqlevels(gr, c("chr2", "CHR1", "chr3")))
    checkIdentical(c("chr2", "CHR1", "chr3"), got)
    got <- seqlevels(suppressWarnings(renameSeqlevels(gr, c(foo="chr2"))))
    checkIdentical(seqlevels(gr), got)
    got <- seqlevels(renameSeqlevels(gr, c(chr2="CHR2")))
    checkIdentical(c("chr1", "CHR2", "chr3"), got)

    ## incomplete rename, order different from seqlevels
    gr <- GRanges(c("chr1", "chr1", "chr2", "chr3"), IRanges(1:4, width=3),
                  seqlengths=c(chr1=10, chr2=20, chr3=30))
    got <- renameSeqlevels(gr, c(chr3="3", chr1="1"))
    checkIdentical(c("1", "1", "chr2", "3"), as.vector(seqnames(got)))
    checkIdentical(structure(c(10L, 20L, 30L), .Names = c("1", "chr2", "3")),
                   seqlengths(got))

    ## expected warning
    txt <- tryCatch(renameSeqlevels(gr, c(chrX="X", chrY="Y", chr1="1")),
                    warning=conditionMessage)
    checkIdentical("invalid seqlevels 'chrX', 'chrY' ignored", txt)
}

test_keepStandardChromosomes <- function()
{      
    gr <- GRanges(c(paste0("chr",1:22), "chrX", "chrY" , "chrM", 
                    "chr1_gl000191_random",
           "chr1_gl000192_random" ,"chr4_ctg9_hap1", "chr4_gl000193_random",
           "chr4_gl000194_random" ,"chr6_apd_hap1"), IRanges(1:31, width=3)) 
    gr2 <- gr
    
    gr <- keepStandardChromosomes(gr)
    checkEquals(25,length(seqlevels(gr)))
           
    gr2 <- keepStandardChromosomes(gr2, species="Homo sapiens")
    checkEquals(25,length(seqlevels(gr2)))
    checkEquals(27,end(gr2[25]))
    
    ## seqlevels common across multiple species. 
    plantgr <- GRanges(c(1:5,"MT","Pltd"), IRanges(1:7,width=5))
    checkEquals(c(1:5,"MT","Pltd"), seqlevels(plantgr))
    
    ## smaller example:
    gr <- GRanges(c("chr1", "chr1", "chr2", "chr3"), IRanges(1:4, width=3))
    gr <- keepStandardChromosomes(gr)
    checkEquals(3, length(seqlevels(gr)))
    checkEquals(c("chr1","chr2","chr3"), seqlevels(gr))
    
    gr <- GRanges(c("chr1", "chr1_gl000192_random", "chrM", "chr1"), 
                  IRanges(1:4, width=3))
    gr <- keepStandardChromosomes(gr) 
    checkEquals(c("chrM","chr1"), seqlevels(gr))
    checkEquals(2, length(seqlevels(gr)))
    
}    

