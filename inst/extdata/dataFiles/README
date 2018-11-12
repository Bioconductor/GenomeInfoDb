The purpose of the files in GenomeInfoDb/inst/extdata/dataFiles is to support
seqlevelsStyle(). seqlevelsStyle() detects a chromosome naming style e.g. UCSC,
NCBI etc.

Many organisms follow the same naming style, whether the style is NCBI, UCSC or
Ensembl.  For example, the chromosome set for Dog (Canis_familiaris.txt) is a
superset of the one for sheep.  seqlevelsStyle() does the right thing on any
vector of chromosome names that follow the same naming style as Dog so it's not
necessary to add a file for Ovis_aries.txt.

Including these base files or supersets eliminates the need to include one
file for each organism we want to support in GenomeInfoDb.
