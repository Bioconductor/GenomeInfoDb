This folder contains tab-delimited files that contain chromosomes info
for a subset of the NCBI assemblies that are registered in the package.

The purpose of these files is to make calls like:

    getChromInfoFromNCBI("GRCh38.p14")

fast and reliable i.e. it should always work, even when offline!

To add a new file, just call:

    saveChromInfoFromNCBI("<assembly-name-or-accession>")

in this folder.

