### =========================================================================
### list_ftp_dir()
### -------------------------------------------------------------------------


### Strip the "protocol://" part (e.g. "ftp://" or "https://") off the
### supplied URL, and add the trailing slash to it if missing.
.normarg_ftp_dir <- function(ftp_dir)
{
    stopifnot(isSingleString(ftp_dir))
    pattern <- "^[[:alpha:]]*://"
    ftp_dir <- sub(pattern, "", ftp_dir)
    ## 'RCurl::getURL(ftp_dir)' will fail if 'ftp_dir' does not have a
    ## trailing slash, so we append the slash but only if it's missing.
    nc <- nchar(ftp_dir)
    last_char <- substr(ftp_dir, nc, nc)
    if (last_char != "/")
        ftp_dir <- paste0(ftp_dir, "/")
    ftp_dir
}

### Wrapper to RCurl::getURL(). Only purpose of this wrapper is to replace
### some uninformative error messages returned by RCurl::getURL() with
### slightly more informative ones.
.getURL2 <- function(url)
{
    stopifnot(isSingleString(url))
    doc <- try(getURL(url), silent=TRUE)
    if (!inherits(doc, "try-error"))
        return(doc)
    condition <- attr(doc, "condition")
    if (inherits(condition, "REMOTE_ACCESS_DENIED")) {
        ## Replace uninformative "Server denied you to change
        ## to the given directory" with more informative message.
        new_msg <- wmsg("Access to \"", url, "\" denied")
        condition$message <- new_msg
    } else if (inherits(condition, "REMOTE_FILE_NOT_FOUND")) {
        ## Replace uninformative "The file does not exist" with
        ## more informative message.
        new_msg <- wmsg("File \"", url, "\" not found on server")
        condition$message <- new_msg
    }
    stop(condition)
}

.make_matrix_data_from_list <- function(x, min_ncol=0L)
{
    stopifnot(is.list(x))
    x_len <- length(x)
    x_lens <- lengths(x)
    ncol <- max(x_lens, min_ncol)
    if (x_len == 0L) {
        ans <- character(0)
    } else {
        y_lens <- ncol - x_lens
        x_seqalong <- seq_along(x)
        f <- rep.int(x_seqalong, y_lens)
        attributes(f) <- list(levels=as.character(x_seqalong), class="factor")
        y <- split(character(length(f)), f)
        collate_subscript <- rep(x_seqalong, each=2L)
        collate_subscript[2L * x_seqalong] <- x_seqalong + x_len
        ans <- unlist(c(x, y)[collate_subscript], use.names=FALSE)
    }
    attr(ans, "ncol") <- ncol
    ans
}

### Return a character matrix with 1 row per entry in 'listing', and at least
### 9 columns (but sometimes more e.g. 11 columns if some entries in 'listing'
### are symlinks). The file names should always be in the 9th column.
.ftp_listing_as_matrix <- function(listing)
{
    split_listing <- strsplit(listing, split="[[:space:]]+")
    ans_data <- .make_matrix_data_from_list(split_listing, min_ncol=9L)
    ans_ncol <- attr(ans_data, "ncol")
    matrix(ans_data, ncol=ans_ncol, byrow=TRUE)
}

### The keys must be URLs to FTP directories that went thru .normarg_ftp_dir(),
### that is, without the "protolol://" part and with a trailing slash.
.cached_ftp_dir_listing <- new.env(parent=emptyenv())

### Uses RCurl::getURL() to retrieve the listing of an FTP dir. The result
### is cached in order to make subsequent calls to list_ftp_dir() on the
### same 'ftp_dir' faster.
### If 'long.listing' is TRUE then return a character matrix with 1 row
### per entry in the FTP dir listing, and at least 9 columns (but sometimes
### more e.g. 11 columns if the listing contains symlinks). The file names
### should always be in the 9th column.
### Note that 'long.listing=TRUE' is similar to option -l of Unix command 'ls'.
list_ftp_dir <- function(ftp_dir, subdirs.only=FALSE, long.listing=FALSE,
                                  recache=FALSE)
{
    ftp_dir <- .normarg_ftp_dir(ftp_dir)
    stopifnot(isTRUEorFALSE(subdirs.only))
    stopifnot(isTRUEorFALSE(long.listing))
    stopifnot(isTRUEorFALSE(recache))

    listing <- .cached_ftp_dir_listing[[ftp_dir]]
    if (is.null(listing) || recache) {
        doc <- .getURL2(ftp_dir)
        listing <- strsplit(doc, "\n", fixed=TRUE)[[1L]]
        .cached_ftp_dir_listing[[ftp_dir]] <- listing
    }

    if (subdirs.only)
        listing <- listing[substr(listing, 1L, 1L) == "d"]
    ans <- .ftp_listing_as_matrix(listing)
    if (!long.listing)
        ans <- ans[ , 9L]
    ans
}

