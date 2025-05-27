This folder contains tab-delimited files that contain chromosomes info
for a subset of the UCSC genomes that are registered in the package.

The purpose of these files is to make calls like:

    getChromInfoFromUCSC("hg38")

fast and reliable i.e. it should always work, even when offline!

To add a new file, just call:

    saveChromInfoFromUCSC("<genome>")

in this folder.

