This folder contains tab-delimited files that contain "assembled molecules
info" for a subset of the UCSC genomes that are registered in the package.

The purpose of these files is to make calls like:

  getChromInfoFromUCSC("hg38", assembled.molecules.only=TRUE)

fast and reliable i.e. it should always work, even when offline!

To add a new file, just call saveAssembledMoleculesInfoFromUCSC() in this
folder.

