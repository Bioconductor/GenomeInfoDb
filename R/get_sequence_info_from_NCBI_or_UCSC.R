### =========================================================================
### get_sequence_info_from_NCBI_or_UCSC()
### -------------------------------------------------------------------------


.get_sequence_info_from_NCBI <- function(genome)
{
    chrom_info <- getChromInfoFromNCBI(genome)
    data.frame(seqnames=chrom_info[ , "SequenceName"],
               seqlengths=chrom_info[ , "SequenceLength"],
               isCircular=chrom_info[ , "circular"],
               genome=genome,
               stringsAsFactors=FALSE)
}

.get_sequence_info_from_UCSC <- function(genome,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    chrom_info <- getChromInfoFromUCSC(genome, goldenPath.url=goldenPath.url)
    data.frame(seqnames=chrom_info[ , "chrom"],
               seqlengths=chrom_info[ , "size"],
               isCircular=chrom_info[ , "circular"],
               genome=genome,
               stringsAsFactors=FALSE)
}

### Return a 4-column data.frame with columns "seqnames" (character),
### "seqlengths" (integer), "isCircular" (logical), and "genome" (character).
get_sequence_info_from_NCBI_or_UCSC <- function(genome,
    goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    if (!isSingleString(genome) || genome == "")
        stop("'genome' must be a single non-empty string")
    NCBI_genomes <- registered_NCBI_genomes()
    UCSC_genomes <- registered_UCSC_genomes()
    if (genome %in% NCBI_genomes[ , "genome"] ||
        genome %in% NCBI_genomes[ , "assembly_accession"])
        return(.get_sequence_info_from_NCBI(genome))
    if (genome %in% UCSC_genomes[ , "genome"])
        return(.get_sequence_info_from_UCSC(genome,
                                            goldenPath.url=goldenPath.url))
    stop(wmsg("\"", genome, "\" is not a registered NCBI or UCSC genome ",
              "(use registered_NCBI_genomes() or registered_UCSC_genomes() ",
              "to list the NCBI or UCSC genome assemblies currently ",
              "registered in the GenomeInfoDb package)"))
}

