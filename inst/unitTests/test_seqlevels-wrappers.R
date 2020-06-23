library(GenomicRanges)

test_keepSeqlevels_GRanges <- function()
{
    gr <- GRanges(c("chr1", "chr4", "chr2", "chr1"), IRanges(1:4, width=3),
                  seqinfo=Seqinfo(c("chr1", "chr2", "chr3", "chr4")))

    ## no pruning
    value <- c("chr4", "chr2", "chr1")
    current <- keepSeqlevels(gr, value)
    checkIdentical(value, seqlevels(current))
    checkIdentical(start(gr), start(current))
    checkIdentical(current, keepSeqlevels(gr, value, pruning.mode="coarse"))
    checkException(keepSeqlevels(gr, value, pruning.mode="fine"), silent=TRUE)
    checkIdentical(current, keepSeqlevels(gr, value, pruning.mode="tidy"))

    ## pruning required
    value <- "chr1"
    checkException(keepSeqlevels(gr, value), silent=TRUE)
    current <- keepSeqlevels(gr, value, pruning.mode="coarse")
    checkIdentical(value, seqlevels(current))
    checkIdentical(c(1L, 4L), start(current))
    checkException(keepSeqlevels(gr, value, pruning.mode="fine"), silent=TRUE)
    checkIdentical(current, keepSeqlevels(gr, value, pruning.mode="tidy"))

    value <- c("chr3", "chr1")
    checkException(keepSeqlevels(gr, value), silent=TRUE)
    current <- keepSeqlevels(gr, value, pruning.mode="coarse")
    checkIdentical(value, seqlevels(current))
    checkIdentical(c(1L, 4L), start(current))
    checkException(keepSeqlevels(gr, value, pruning.mode="fine"), silent=TRUE)
    checkIdentical(current, keepSeqlevels(gr, value, pruning.mode="tidy"))

    value <- character(0)
    checkException(keepSeqlevels(gr, value), silent=TRUE)
    current <- keepSeqlevels(gr, value, pruning.mode="coarse")
    checkIdentical(value, seqlevels(current))
    checkIdentical(integer(0), start(current))
    checkException(keepSeqlevels(gr, value, pruning.mode="fine"), silent=TRUE)
    checkIdentical(current, keepSeqlevels(gr, value, pruning.mode="tidy"))

    ## invlide seqlevels
    value <- c("chr3", "chrX", "chr1", "chrY")
    checkException(keepSeqlevels(gr, value), silent=TRUE)
    checkException(keepSeqlevels(gr, value, pruning.mode="coarse"), silent=TRUE)
    checkException(keepSeqlevels(gr, value, pruning.mode="fine"), silent=TRUE)
    checkException(keepSeqlevels(gr, value, pruning.mode="tidy"), silent=TRUE)
}

test_dropSeqlevels_GRanges <- function()
{
    gr <- GRanges(c("chr1", "chr4", "chr2", "chr1"), IRanges(1:4, width=3),
                  seqinfo=Seqinfo(c("chr1", "chr2", "chr3", "chr4")))

    ## no pruning
    value <- c("chrX", "chr3", "chr3")  # duplicates are ok
    current <- dropSeqlevels(gr, value)
    checkIdentical(c("chr1", "chr2", "chr4"), seqlevels(current))
    checkIdentical(start(gr), start(current))
    checkIdentical(current, dropSeqlevels(gr, value, pruning.mode="coarse"))
    checkException(dropSeqlevels(gr, value, pruning.mode="fine"), silent=TRUE)
    checkIdentical(current, dropSeqlevels(gr, value, pruning.mode="tidy"))

    ## pruning required
    value <- c("chrX", "chr3", "chr1", "chr1")  # duplicates are ok
    checkException(dropSeqlevels(gr, value), silent=TRUE)
    current <- dropSeqlevels(gr, value, pruning.mode="coarse")
    checkIdentical(c("chr2", "chr4"), seqlevels(current))
    checkIdentical(2:3, start(current))
    checkException(dropSeqlevels(gr, value, pruning.mode="fine"), silent=TRUE)
    checkIdentical(current, dropSeqlevels(gr, value, pruning.mode="tidy"))

    value <- character(0)
    checkIdentical(gr, dropSeqlevels(gr, value))
    checkIdentical(gr, dropSeqlevels(gr, value, pruning.mode="coarse"))
    checkException(dropSeqlevels(gr, value, pruning.mode="fine"), silent=TRUE)
    checkIdentical(gr, dropSeqlevels(gr, value, pruning.mode="tidy"))

    value <- c("chrX", rev(seqlevels(gr)))
    checkException(dropSeqlevels(gr, value), silent=TRUE)
    current <- dropSeqlevels(gr, value, pruning.mode="coarse")
    checkIdentical(character(0), seqlevels(current))
    checkIdentical(integer(0), start(current))
    checkException(dropSeqlevels(gr, value, pruning.mode="fine"), silent=TRUE)
    checkIdentical(current, dropSeqlevels(gr, value, pruning.mode="tidy"))
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
    
    gr <- keepStandardChromosomes(gr, pruning.mode="coarse")
    checkEquals(25,length(seqlevels(gr)))
    
    ## want to test if sorted seqlevels are being returned.
    checkEquals(seqlevels(gr), c(paste0("chr",1:22), "chrX", "chrY" , "chrM"))
           
    gr2 <- keepStandardChromosomes(gr2, species="Homo sapiens",
                                        pruning.mode="coarse")
    checkEquals(25,length(seqlevels(gr2)))
    checkEquals(27,end(gr2[25]))
    
    gr <- GRanges(c("chr1", "chr2", "chr3L", "3L"), 
                  IRanges(1:4, width=3))
    gr <- keepStandardChromosomes(gr, species="Homo sapiens",
                                      pruning.mode="coarse") 
    checkEquals(2, length(seqlevels(gr)))
   
    ## drop scaffolds - eg1. 
    gr <- GRanges(c("chr1", "chr1_gl000192_random", "chrM", "chr1"), 
                  IRanges(1:4, width=3))
    gr <- keepStandardChromosomes(gr, pruning.mode="coarse")
    checkEquals(c("chr1","chrM"), seqlevels(gr))
    checkEquals(2, length(seqlevels(gr)))
    
    ## seqlevels not supported by GenomeInfodb.  
    gr <- GRanges("chr4_ctg9_hap1", IRanges(1, 5))
    checkEquals(
        length(seqlevels(keepStandardChromosomes(gr, pruning.mode="coarse"))),
        0
    )
    checkException(seqlevelsStyle(gr), silent=TRUE)
    
    ## drop seqlevels not supported by GenomeInfoDb
    plantgr <- GRanges(c(41:45,"MT","Pltd"), IRanges(1:7,width=5))
    plantgr <- keepStandardChromosomes(plantgr, pruning.mode="coarse")
    checkEquals(c("MT","Pltd"), seqlevels(plantgr))
    
    ## no seqlevels in object
    checkEquals(0,length(seqlevels(keepStandardChromosomes(GRanges()))))
    
    ## mixed seqlevels; no species
    gr <- GRanges(c("chr8", "chr2", "foo", "chr3L"), IRanges(1:4, width=3))
    ans <- keepStandardChromosomes(gr, pruning.mode="coarse")
    checkEquals(seqlevels(ans), c("chr8", "chr2", "chr3L"))
    gr <- GRanges(c("chr8", "chr2", "chr3R", "chr3L"), IRanges(1:4, width=3))
    ans <- keepStandardChromosomes(gr)
    checkEquals(seqlevels(ans), c("chr8", "chr2", "chr3R", "chr3L"))
    ## mixed seqlevels; species
    gr <- GRanges(c("chr8", "chr2", "foo", "chr3L"), IRanges(1:4, width=3))
    ans <- keepStandardChromosomes(gr, "Homo sapiens", pruning.mode="coarse")
    checkEquals(seqlevels(ans), c("chr8", "chr2"))
    fly <- "Drosophila melanogaster"
    ans <- keepStandardChromosomes(gr, fly, pruning.mode="coarse")
    checkEquals(seqlevels(ans), "chr3L")

    # match multiple styles
    gr <- GRanges( c("1","1","X","FOO"), IRanges(start=1:4,width=2))
    ans <- keepStandardChromosomes(gr, pruning.mode="coarse")
    checkEquals(seqlevels(ans), c("1", "X"))
    ans <- keepStandardChromosomes(gr, species="Homo sapiens",
                                       pruning.mode="coarse")
    checkEquals(seqlevels(ans), c("1", "X"))
}    

