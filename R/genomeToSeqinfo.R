### =========================================================================
### genomeToSeqinfo()
### -------------------------------------------------------------------------


.UCSC_genome_to_Seqinfo <- function(genome, circ_seqs,
    goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
{
    chromInfo_table <- fetch_ChromInfo_from_UCSC(genome, goldenPath_url)
    ans_seqnames <- as.character(chromInfo_table[ , "chrom"])
    oo <- order(rankSeqlevels(ans_seqnames))
    ans_seqnames <- ans_seqnames[oo]
    ans_seqlengths <- as.integer(chromInfo_table[ , "size"])[oo]
    ans_isCircular <- logical(nrow(chromInfo_table))
    if (!is.null(circ_seqs)) {
        circ_idx <- match(circ_seqs, ans_seqnames)
        if (any(is.na(circ_idx)))
            stop("'circ_seqs' contains sequence names not in ",
                 genome, " genome")
        ans_isCircular[circ_idx] <- TRUE
    }
    Seqinfo(ans_seqnames, ans_seqlengths, ans_isCircular, genome)
}

### Only supports UCSC genomes for now (the same genomes that are supported
### by fetchExtendedChromInfoFromUCSC()).
### NOT exported.
genomeToSeqinfo <- function(genome)
{
    if (!isSingleString(genome) || genome == "")
        stop("'genome' must be a single non-empty string")
    idx <- match(genome, names(SUPPORTED_UCSC_GENOMES))
    if (is.na(idx))
        stop("genome \"", genome, "\" is not supported")
    circ_seqs <- SUPPORTED_UCSC_GENOMES[[idx]]$circular
    .UCSC_genome_to_Seqinfo(genome, circ_seqs)
}

