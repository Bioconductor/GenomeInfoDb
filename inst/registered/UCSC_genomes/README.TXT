This folder contains "registration files" for UCSC genomes.

These files are used by the getChromInfoFromUCSC() function defined in
the GenomeInfoDb package.

There must be one file per genome.

Each file must be an R script (.R extension) and its name must be the name
of the genome (e.g. 'danRer11.R').

The script should be able to work as a standalone script so should
explicitly load packages if needed (e.g. with 'library(IRanges)'
and/or 'library(GenomeInfoDb)').

At a minimum, the script must define the 4 following variables:

  o GENOME:              Single non-empty string.

  o ORGANISM:            Single non-empty string.

  o ASSEMBLED_MOLECULES: Character vector with no NAs, no empty strings,
                         and no duplicates.

  o CIRC_SEQS:           Character vector (subset of ASSEMBLED_MOLECULES).

Additionally, it can also define any of the following variables:

  o FETCH_ORDERED_CHROM_SIZES: Function with 1 argument. Must return a 2-column
                         data frame with columns "chrom" (character)
                         and "size" (integer). Rows must be in "canonical"
                         chromosome order. A requirement is that the assembled
                         molecules come first.
                         Note that defining this function is not needed if
                         all the sequences in the genome are assembled
                         molecules. See for example registration files for
                         Worm (ce*.R files).

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

The script can define its own top-level variables and functions, but, by
convention, their names should be in lower case and start with a dot.

See the files in this folder for numerous examples.

Here is how to perform some basic testing of a new registration file:

  1. In a **fresh** R session, use source() to source the new file.
     This has the effect of executing the code in the script (alternatively
     you can copy-paste the content of the script in your session).
     Note that a registration file is required to be a **standalone** R
     script. This means that we should be able to source it in a fresh R
     session, and it should just work (granted that all the required packages
     are installed). "Just work" here means that all the top-level variables
     defined in the script will get defined in your session (you should see
     them with ls()).

  2. Check the values of the 4 mandatory variables: GENOME, ORGANISM,
     ASSEMBLED_MOLECULES, and CIRC_SEQS.

  3. Call FETCH_ORDERED_CHROM_SIZES() and make sure it behaves has expected.
     Does the returned data frame has its rows in the expected order?

  4. Install GenomeInfoDb (with the new registration file in it), start R,
     load the package, and try to call registered_UCSC_genomes(). The
     returned data frame should now have an entry for the new genome.
     Make sure that all the fields in the new entry look as expected.
     Then call getChromInfoFromUCSC() on the new genome. Do the "assembled"
     and "circular" columns look as expected? Try with and without
     setting the 'assembled.molecules.only' argument to TRUE.

  5. If you've defined NCBI_LINKER in the registration file: try to call
     getChromInfoFromUCSC() with 'map.NCBI' set to TRUE.

  6. If you've defined ENSEMBL_LINKER in the registration file: try to call
     getChromInfoFromUCSC() with 'add.ensembl.col' set to TRUE.

