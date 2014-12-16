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
                 "\"", genome, "\" genome")
        ans_isCircular[circ_idx] <- TRUE
    }
    Seqinfo(ans_seqnames, ans_seqlengths, ans_isCircular, genome)
}

#.NCBI_genome_to_Seqinfo <- function(refseq_assembly_id, AssemblyUnits,
#                                                        circ_seqs)
#{
#    assembly_report <- fetch_assembly_report(refseq_assembly_id,
#                                             AssemblyUnits=AssemblyUnits)
#    ans_seqnames <- as.character(assembly_report[ , "SequenceName"])
#
#}
#
#SUPPORTED_NCBI_GENOMES <- list(
#    GRCh38=
#        list(refseq_assembly_id="GCF_000001405.26", circular="MT")
#)

### Only supports UCSC genomes for now (the same genomes that are supported
### by fetchExtendedChromInfoFromUCSC()).
### NOT exported.
genomeToSeqinfo <- function(genome)
{
    if (!isSingleString(genome) || genome == "")
        stop("'genome' must be a single non-empty string")
    idx <- match(genome, names(SUPPORTED_UCSC_GENOMES))
    if (!is.na(idx)) {
        supported_genome <- SUPPORTED_UCSC_GENOMES[[idx]]
        circ_seqs <- supported_genome$circular
        return(.UCSC_genome_to_Seqinfo(genome, circ_seqs))
    }
    #idx <- match(genome, names(SUPPORTED_NCBI_GENOMES))
    #if (!is.na(idx)) {
    #    supported_genome <- SUPPORTED_NCBI_GENOMES[[idx]]
    #    refseq_assembly_id <- supported_genome$refseq_assembly_id
    #    AssemblyUnits <- supported_genome$AssemblyUnits
    #    circ_seqs <- supported_genome$circular
    #    return(.NCBI_genome_to_Seqinfo(refseq_assembly_id, AssemblyUnits,
    #                                                       circ_seqs,))
    #}
    stop("genome \"", genome, "\" is not supported")
}

