require("GenomeInfoDb") || stop("unable to load GenomeInfoDb package")
require("RUnit") || stop("unable to load RUnit package")
GenomeInfoDb:::.test()
