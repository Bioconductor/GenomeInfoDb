### =========================================================================
### Miscellaneous low-level utils
### -------------------------------------------------------------------------


### Global character vector to hold default names for circular sequences.
### This is exported!
DEFAULT_CIRC_SEQS <- c(
    ## Mitochondrial genome
    "chrM", "MT", "MtDNA", "mit", "Mito", "mitochondrion",
    "dmel_mitochondrion_genome",
    ## Chloroplast genome
    "Pltd", "ChrC", "Pt", "chloroplast", "Chloro",
    ## Plasmid (yeast)
    "2micron", "2-micron", "2uM"
)

### AFAIK UCSC doesn't flag circular sequences.
### As of Sep 21, 2010 (Ensembl release 59), Ensembl was still not flagging
### circular sequences in their db (see this thread for the details
### http://lists.ensembl.org/pipermail/dev/2010-September/000139.html),
### NOT exported but used in the GenomicFeatures package.
make_circ_flags_from_circ_seqs <- function(seqlevels, circ_seqs=NULL)
{
    if (!is.character(seqlevels))
        stop(wmsg("'seqlevels' must be a character vector"))
    if (is.null(circ_seqs)) {
        ## The user did NOT specify the 'circ_seqs' argument.
        seqlevels <- tolower(seqlevels)
        circ_seqs <- tolower(DEFAULT_CIRC_SEQS)
        circ_flags <- rep.int(NA, length(seqlevels))
        circ_flags[seqlevels %in% circ_seqs] <- TRUE
    } else {
        ## The user specified the 'circ_seqs' argument.
        if (!is.character(circ_seqs) ||
            any(circ_seqs %in% c(NA_character_, "")))
            stop(wmsg("'circ_seqs' must be a character vector with no NAs ",
                      "and no empty strings"))
        bad_circ_seqs <- setdiff(circ_seqs, seqlevels)
        if (length(bad_circ_seqs) != 0L) {
            in1string <- paste0(bad_circ_seqs, collapse=", ")
            stop(wmsg("'circ_seqs' contains unrecognized chromosome names: ",
                      in1string))
        }
        circ_flags <- seqlevels %in% circ_seqs
    }
    circ_flags
}

### Note that, strictly speaking, mergeNamedAtomicVectors() is not
### commutative, i.e., in general 'z1 <- mergeNamedAtomicVectors(x, y)' is
### not identical to 'z2 <- mergeNamedAtomicVectors(y, x)'. However 'z1' and
### 'z2' are both guaranteed to have unique names and to contain the same set
### of name/value pairs (but typically not in the same order).
mergeNamedAtomicVectors <- function(x, y, what=c("key", "values"))
{
    if (!is.atomic(x) || !is.atomic(y) || typeof(x) != typeof(y))
        stop("'x' and 'y' must be atomic vectors of the same type")
    x_names <- names(x)
    y_names <- names(y)
    if (is.null(x_names) || is.null(y_names))
        stop("'x' and 'y' must have names")
    ans_names <- union(x_names, y_names)
    if (any(ans_names %in% c(NA_character_, "")))
        stop("some names in 'x' or 'y' are NA or the empty string")
    ## Note sure why but subsetting by name is *very* slow when the character
    ## vector used as subscript contains a lot of "invalid" names. I already
    ## reported this issue twice on the R-devel mailing list (in July 2010 and
    ## May 2013) with no answer so far.
    #ans <- x[ans_names]  # very slow :-(
    ans <- x[match(ans_names, x_names)]  # much faster :-)
    ## Some of 'ans' names can be NA. This is because subsetting a named
    ## vector 'x' by a character vector 'i' returns a vector whose names
    ## are 'i' except for the values in 'i' that are not in 'names(x)'.
    ## So we need to fix this.
    names(ans) <- ans_names
    #ans2 <- y[ans_names]  # very slow :-(
    ans2 <- y[match(ans_names, y_names)]  # much faster :-)
    idx <- which(ans != ans2)
    if (length(idx) != 0L) {
        msg <- c(what[1L], ifelse(length(idx) >= 2, "s", ""), " ",
                 paste(ans_names[idx], collapse=", "), " ",
                 ifelse(length(idx) >= 2, "have", "has"),
                 " incompatible ", what[2L], ":\n  - in 'x': ",
                 paste(ans[idx], collapse=", "), "\n  - in 'y': ",
                 paste(ans2[idx], collapse=", "))
        stop(msg)
    }
    idx <- is.na(ans) & !is.na(ans2)
    ans[idx] <- ans2[idx]
    ans
}

### Uses RCurl to access and list the content of an FTP dir.
list_ftp_dir <- function(url)
{
    doc <- getURL(url)  # from RCurl package
    listing <- strsplit(doc, "\n", fixed=TRUE)[[1L]]
    ## Keep field no. 8 only
    pattern <- paste(c("^", rep.int("[^[:space:]]+[[:space:]]+", 8L)),
                     collapse="")
    listing <- sub(pattern, "", listing)
    sub("[[:space:]].*$", "", listing)
}

