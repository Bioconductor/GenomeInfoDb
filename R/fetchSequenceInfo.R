### =========================================================================
### fetchSequenceInfo()
### -------------------------------------------------------------------------


.fetch_sequence_info_for_UCSC_genome <- function(genome,
    goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
{
    ext_chrominfo <- fetchExtendedChromInfoFromUCSC(genome,
                         goldenPath_url=goldenPath_url, quiet=TRUE)
    data.frame(seqnames=ext_chrominfo[ , "UCSC_seqlevels"],
               seqlengths=ext_chrominfo[ , "UCSC_seqlengths"],
               is_circular=ext_chrominfo[ , "circular"],
               genome=genome,
               stringsAsFactors=FALSE)
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
    if (!is.na(idx))
        return(.fetch_sequence_info_for_UCSC_genome(genome))
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

