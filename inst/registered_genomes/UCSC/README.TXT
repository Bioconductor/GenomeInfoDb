This folder contains "registration files" for UCSC genome assemblies.

These files are used by the getChromInfoFromUCSC() function defined in
the GenomeInfoDb package.

There must be one file per assembly.

Each file must be an R script (.R extension) and its name must be the name
of the genome assembly (e.g. 'danRer11.R').

The script should be able to work as a standalone script so should
explicitly load packages if needed (e.g. with 'library(IRanges)'
and/or 'library(GenomeInfoDb)').

At a minimum, the script must define the 4 following variables:

  o GENOME:              Single non-empty string.

  o ORGANISM:            Single non-empty string.

  o ASSEMBLED_MOLECULES: Character vector with no NAs, no empty strings,
                         and no duplicates.

  o CIRC_SEQS:           Character vector (subset of ASSEMBLED_MOLECULES).

Additionally, it can also define any (or both) of the following variables:

  o GET_CHROM_SIZES:     Function with 1 argument. Must return a 2-column
                         data frame with columns "chrom" (character)
                         and "size" (integer).

  o NCBI_LINKER:         Named list.

    Valid NCBI_LINKER components:
    - assembly_accession: single non-empty string.
    - AssemblyUnits: character vector.
    - special_mappings: named character vector.
    - unmapped_seqs: named list of character vectors.
    - drop_unmapped: TRUE or FALSE.

  o ENSEMBL_LINKER:      Single non-empty string (can only be "ucscToEnsembl"
                         or "chromAlias" at the moment).

All the above variables are recognized by getChromInfoFromUCSC(). They
must be defined at the top-level of the script and their names must be
in UPPER CASE.

The script can define its own top-level variables and functions but by
convention their name should be in lower case and start with a dot.

See the files in this folder for numerous examples.

