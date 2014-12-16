### =========================================================================
### fetchSequenceInfo()
### -------------------------------------------------------------------------


.fetch_sequence_info_for_UCSC_genome <- function(genome, circ_seqs,
    goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
{
    chromInfo_table <- fetch_ChromInfo_from_UCSC(genome, goldenPath_url)
    ans_seqnames <- as.character(chromInfo_table[ , "chrom"])
    ans_seqlengths <- as.integer(chromInfo_table[ , "size"])
    ans_is_circular <- logical(nrow(chromInfo_table))
    if (!is.null(circ_seqs)) {
        circ_idx <- match(circ_seqs, ans_seqnames)
        if (any(is.na(circ_idx)))
            stop("'circ_seqs' contains sequence names not in ",
                 "\"", genome, "\" genome")
        ans_is_circular[circ_idx] <- TRUE
    }
    sequence_info <- data.frame(seqnames=ans_seqnames,
                                seqlengths=ans_seqlengths,
                                is_circular=ans_is_circular,
                                genome=genome,
                                stringsAsFactors=FALSE)
    sequence_info[order(rankSeqlevels(ans_seqnames)), , drop=FALSE]
}

#.fetch_sequence_info_for_NCBI_genome <- function(refseq_assembly_id,
#                                                 AssemblyUnits,
#                                                 circ_seqs)
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

### Returns a data frame.
### Only supports UCSC genomes for now (the same genomes that are supported
### by fetchExtendedChromInfoFromUCSC()).
### NOT exported.
fetchSequenceInfo <- function(genome)
{
    if (!isSingleString(genome) || genome == "")
        stop("'genome' must be a single non-empty string")
    idx <- match(genome, names(SUPPORTED_UCSC_GENOMES))
    if (!is.na(idx)) {
        supported_genome <- SUPPORTED_UCSC_GENOMES[[idx]]
        circ_seqs <- supported_genome$circular
        return(.fetch_sequence_info_for_UCSC_genome(genome, circ_seqs))
    }
    #idx <- match(genome, names(SUPPORTED_NCBI_GENOMES))
    #if (!is.na(idx)) {
    #    supported_genome <- SUPPORTED_NCBI_GENOMES[[idx]]
    #    refseq_assembly_id <- supported_genome$refseq_assembly_id
    #    AssemblyUnits <- supported_genome$AssemblyUnits
    #    circ_seqs <- supported_genome$circular
    #    return(.fetch_sequence_info_for_NCBI_genome(refseq_assembly_id,
    #                                                AssemblyUnits,
    #                                                circ_seqs))
    #}
    stop("genome \"", genome, "\" is not supported")
}

