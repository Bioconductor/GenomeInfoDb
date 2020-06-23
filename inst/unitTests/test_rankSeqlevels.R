test_orderSeqlevels <- function()
{
    checkIdentical(integer(), orderSeqlevels())
    checkException(orderSeqlevels(c(1,2,3,4,5)))

    seqnames <- c("chr1","chr9", "chr2", "chr3", "chr10")
    checkEquals(c(1L, 3L, 4L, 2L, 5L), orderSeqlevels(seqnames) )
}

