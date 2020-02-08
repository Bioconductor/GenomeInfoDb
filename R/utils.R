### =========================================================================
### Miscellaneous low-level utils
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### 'x' must be a matrix-like or data-frame-like object.
drop_cols <- function(x, colnames)
{
    stopifnot(is.character(colnames), length(colnames) != 0L)
    drop_idx <- match(colnames, colnames(x))
    stopifnot(!anyNA(drop_idx))
    x[ , -drop_idx, drop=FALSE]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### JOIN two data.frames
###
### In typical base R fashion, merge.data.frame() does all kinds of weird
### things with the column names and is very inefficient. Plus, independently
### of how 'sort' is set, the order of the rows in the result does not make
### sense. The order of the columns doesn't either.
###
### join_dfs() below does NOT mess around with the column names, preserves
### the order of the rows of the left data.frame, does not shuffle the
### columns around (left and right columns stay on the left and right,
### respectively, and in their original order), and is about 3x-4x faster
### than merge.data.frame()!

.do_join <- function(Ldf, Rdf, L2R)
{
    #stopifnot(is.integer(L2R),
    #          length(L2R) == nrow(Ldf),
    #          all(L2R >= 1L, na.rm=TRUE),
    #          all(L2R <= nrow(Rdf), na.rm=TRUE))
    cbind(Ldf, S4Vectors:::extract_data_frame_rows(Rdf, L2R))
}

### Performs an SQL INNER JOIN by default. Also by default the right column
### of the join ('Rcolname') is not included in the result.
### IMPORTANT NOTE: It will mimic the behavior of an SQL JOIN (INNER or LEFT)
### **only** if the right column contains unique values! (UNIQUE constraint
### in SQL). If not, then values in the left column ('Lcolname') will only
### get mapped to their first match in the right column ('Rcolname').
### Comparison with merge.data.frame():
### - INNER JOIN (join_dfs is about 3x faster):
###     join_dfs(df1, df2, Lcolname, Rcolname)
###     merge(df1, df2, by.x=Lcolname, by.y=Rcolname, sort=FALSE)
### - LEFT JOIN (join_dfs is about 4x faster):
###     join_dfs(df1, df2, Lcolname, Rcolname, left.join=TRUE)
###     merge(df1, df2, by.x=Lcolname, by.y=Rcolname, sort=FALSE, all.x=TRUE)
join_dfs <- function(Ldf, Rdf, Lcolname, Rcolname,
                     left.join=FALSE, keep.Rcol=FALSE)
{
    stopifnot(is.data.frame(Ldf),
              is.data.frame(Rdf),
              isSingleString(Lcolname),
              isSingleString(Rcolname))
    ## Values in the left column only get mapped to their first match in
    ## the right column.
    L2R <- match(Ldf[ , Lcolname], Rdf[ , Rcolname])
    if (!left.join) {
        ## Drop rows in 'Ldf' that are not mapped.
        drop_idx <- which(is.na(L2R))
        if (length(drop_idx) != 0L) {
            Ldf <- S4Vectors:::extract_data_frame_rows(Ldf, -drop_idx)
            L2R <- L2R[-drop_idx]
        }
    }
    if (!keep.Rcol)
        Rdf <- drop_cols(Rdf, Rcolname)
    .do_join(Ldf, Rdf, L2R)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### solid_match() and solid_match2()
###

.one_to_many_msg <- function(culprits, from_what, to_what)
{
    culprits_names <- names(culprits)
    if (!is.null(culprits_names))
        culprits <- culprits_names
    in1string <- paste0(unique(culprits), collapse=", ")
    c(from_what, "(s) matched to more than 1 ", to_what, ": ", in1string)
}

### Like base::match() but raises an error if some elements in 'x' can be
### matched to more than one element in 'table' (one-to-many mapping).
### IMPORTANT NOTE: The fast implementation below doesn't work if 'x'
### contains duplicates so we check this. This means that we cannot use
### it in join_dfs() above!
solid_match <- function(x, table,
                        x_what="'x' element", table_what="'table' element")
{
    ## We only support atomic vectors (this is unlike base::match()
    ## where 'x' and 'table' each can be a NULL or a list).
    stopifnot(is.vector(x), is.atomic(x), !anyDuplicated(x),
              is.vector(table), is.atomic(table))
    revm <- match(table, x)  # reverse match
    ambig_idx <- which(duplicated(revm, incomparables=NA_integer_))
    if (length(ambig_idx) != 0L) {
        uidx <- unique(revm[ambig_idx])
        stop(wmsg(.one_to_many_msg(x[uidx], x_what, table_what)))
    }
    ans <- rep.int(NA_integer_, length(x))
    ok <- !is.na(revm)
    ans[revm[ok]] <- which(ok)
    ans
}

### Like base::match() but raises an error if the direct or reverse mapping
### is one-to-many.
solid_match2 <- function(x, table,
                         x_what="'x' element", table_what="'table' element")
{
    ## We only support atomic vectors (this is unlike base::match()
    ## where 'x' and 'table' each can be a NULL or a list).
    stopifnot(is.vector(x), is.atomic(x),
              is.vector(table), is.atomic(table))
    hits <- findMatches(x, table)
    q_hits <- queryHits(hits)
    s_hits <- subjectHits(hits)
    q_ambig_idx <- which(duplicated(q_hits))
    if (length(q_ambig_idx) != 0L) {
        uidx <- unique(q_hits[q_ambig_idx])
        stop(wmsg(.one_to_many_msg(x[uidx], x_what, table_what)))
    }
    s_ambig_idx <- which(duplicated(s_hits))
    if (length(s_ambig_idx) != 0L) {
        uidx <- unique(s_hits[s_ambig_idx])
        stop(wmsg(.one_to_many_msg(table[uidx], table_what, x_what)))
    }
    ans <- rep.int(NA_integer_, length(x))
    ans[q_hits] <- s_hits
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other stuff
###

### Uses RCurl to access and list the content of an FTP dir.
list_ftp_dir <- function(url, subdirs.only=FALSE)
{
    doc <- getURL(url)
    listing <- strsplit(doc, "\n", fixed=TRUE)[[1L]]
    if (subdirs.only)
        listing <- listing[substr(listing, 1L, 1L) == "d"]
    ## Keep field no. 8 only
    pattern <- paste(c("^", rep.int("[^[:space:]]+[[:space:]]+", 8L)),
                     collapse="")
    listing <- sub(pattern, "", listing)
    sub("[[:space:]].*$", "", listing)
}

### Provides a simpler interface to read.table().
.simple_read_table <- function(file,
                               header=FALSE, colnames=NULL, col2class=NULL,
                               nrows=-1L, skip=0L,
                               sep="\t", quote="", comment.char="")
{
    if (is.null(col2class))
        col2class <- NA
    ## Prepare args to pass to read.table().
    args <- list(file, header=header,
                 sep=sep, quote=quote, na.strings=c("NA", "na"),
                 colClasses=col2class, nrows=nrows, skip=skip,
                 comment.char=comment.char, stringsAsFactors=FALSE)
    if (!is.null(colnames))
        args$col.names <- colnames
    do.call(read.table, args)
}

### Same interface as .simple_read_table() above.
fetch_table_from_url <- function(url, ...)
{
    destfile <- tempfile()
    download.file(url, destfile, quiet=TRUE)
    ans <- .simple_read_table(destfile, ...)
    unlink(destfile)
    ans
}

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

