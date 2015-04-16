### =========================================================================
### rankSeqlevels()
### -------------------------------------------------------------------------
###
### Assign a unique ID to each unique sequence name passed in 'seqnames'.
### The returned IDs are guaranteed to span 1:N where N is the number of
### unique sequence names in 'seqnames'.
### Also the function tries hard to assign IDs in a way that is consistent
### with a "good looking" order of the sequence names. This "good looking"
### order is roughly defined by the following (complicated and arbitrary)
### set of rules (rules apply in the order shown below):
###
###   1. Every name should fall into exactly 1 of the 5 following "super
###      groups":
###        (A) starts with CHR
###        (B) starts with chr
###        (C) starts with CH
###        (D) starts with ch
###        (E) anything else
###      Names in early super groups are ranked before names in late super
###      groups.
###
###   2. Within each super group, and after the prefix corresponding to the
###      super group has been dropped (nothing is dropped for super group (E)),
###      every name should fall into exactly 1 of the 18 following groups:
###        (a) roman number
###        (b) "short" arabic integer number (i.e. 6 digits or less) possibly
###            with A, a, B, b, L, or R suffix
###        (c) W
###        (d) Z
###        (e) X
###        (f) Y
###        (g) U
###        (h) M
###        (i) MT
###        (j) "short" arabic integer number (i.e. 6 digits or less) "followed
###            by something" (not A, a, B, b, L, or R)
###        (k) W "followed by something"
###        (l) Z "followed by something"
###        (m) X "followed by something"
###        (n) Y "followed by something"
###        (o) U "followed by something"
###        (p) M "followed by something"
###        (q) MT "followed by something"
###        (r) anything else
###      Names in early groups are ranked before names in late groups.
###
###   3. A name in group (b) with A, a, B, b, L, or R suffix, is ranked
###      right after the name obtained by dropping the suffix.
###
###   4. In groups (k-r), ties are broken by looking at the "followed by
###      something" part (or at the entire name for group (r)): collation
###      defined by LC_COLLATE set to C applies.
###
### 'X.is.sexchrom' lets the user control whether X refers to the sexual
### chromosome or to chromosome with roman number X.
###
### Yes, an ugly and messy function, sorry...
###
### NOTE: rankSeqlevels() was successfully tested on the BSgenome data
### packages for hg19, mm10, ce2, dm3, sacCer1, sacCer2, sacCer3 and rheMac2
### i.e. the IDs returned on the seqnames defined in those packages match the
### ranks of the seqnames.
### For example, for hg19, 'rankSeqlevels(seqlevels(Hsapiens))' is identical
### to 'seq_along(seqlevels(Hsapiens)))'.
### TODO: Add unit test for rankSeqlevels().


### Some simple helpers for low-level string manipulation.

.hasPrefix <- function(x, prefix)
    substr(x, start=1L, stop=nchar(prefix)) == prefix

.dropPrefix <- function(x, nchar)
    substr(x, start=nchar+1L, stop=nchar(x))

isRoman <- function(x)
{
    suppressWarnings(roman <- utils:::.roman2numeric(x))
    ans <- logical(length(x))
    ans[!is.na(roman) & toupper(x) == x] <- TRUE
    ans
}
    
.REGEXP0 <- "[1-9][0-9]*"
.REGEXP1 <- "[0-9]*"

.getNbPart <- function(x)
{
    pattern <- paste0("^(", .REGEXP1, ")(.*)$")
    sub(pattern, "\\1", x)
}

.getPostNbPart <- function(x)
{
    pattern <- paste0("^(", .REGEXP1, ")(.*)$")
    sub(pattern, "\\2", x)
}

.isShortNb <- function(x, abc="")
{
    pattern <- paste0("^", .REGEXP0, abc, "$")
    is_nb <- grepl(pattern, x)
    prefix <- .getNbPart(x)
    is_nb & nchar(prefix) <= 6L
}

### rankSeqlevels()

rankSeqlevels <- function(seqnames, X.is.sexchrom=NA)
{
    if (is.character(seqnames))
        seqnames <- factor(seqnames)
    else if (!is.factor(seqnames))
        stop("'seqnames' must be a character vector or factor")
    if (!is.logical(X.is.sexchrom) || length(X.is.sexchrom) != 1L)
        stop("'X.is.sexchrom' must be a single logical")
    seqlevels <- levels(factor(seqnames))  # unique seqnames
    
    ## Set LC_COLLATE to C so our calls to as.integer(factor(...)) below
    ## return the same thing for everybody (i.e. for every user on any machine
    ## in any country).
    prev_locale <- Sys.getlocale("LC_COLLATE")
    Sys.setlocale("LC_COLLATE", "C")
    on.exit(Sys.setlocale("LC_COLLATE", prev_locale))
    
    ## Provisional ids.
    prov_ids <- rep.int(NA_integer_, length(seqlevels))
    last_prov_id <- 0L
    
    ## 'i' indices of elements in 'prov_ids' to set.
    ## 'ints' integer vector of length > 0 with no NAs. Is recycled to the
    ## length of 'i'.
    makeAndAssignProvIds <- function(i, ints=0L)
    {
        new_prov_ids <- last_prov_id + ints - min(ints) + 1L
        prov_ids2 <- prov_ids
        prov_ids2[i] <- new_prov_ids
        last_prov_id2 <- max(new_prov_ids)
        assign("prov_ids", prov_ids2, inherits=TRUE)
        assign("last_prov_id", last_prov_id2, inherits=TRUE)
    }
    
    assignProvIdsForSuperGroup <- function(seqlevels, prefix)
    {
        sgidx <- which(is.na(prov_ids) & .hasPrefix(seqlevels, prefix))
        if (length(sgidx) == 0L)
            return()
        sgsuffix <- seqlevels[sgidx]
        sgsuffix <- .dropPrefix(sgsuffix, nchar(prefix))
        is_nb <- .isShortNb(sgsuffix)
        is_nbA <- .isShortNb(sgsuffix, abc="A")
        is_nba <- .isShortNb(sgsuffix, abc="a")
        is_nbB <- .isShortNb(sgsuffix, abc="B")
        is_nbb <- .isShortNb(sgsuffix, abc="b")
        is_nbL <- .isShortNb(sgsuffix, abc="L")
        is_nbR <- .isShortNb(sgsuffix, abc="R")
        is_nb_with_suffix <- is_nb | is_nbA | is_nba | is_nbB | is_nbb |
                                                       is_nbL | is_nbR
        is_nbxxx <- .isShortNb(sgsuffix, abc=".*") & !is_nb_with_suffix
        is_W <- sgsuffix == "W"
        is_Z <- sgsuffix == "Z"
        is_X <- sgsuffix == "X"
        is_Y <- sgsuffix == "Y"
        if (is.na(X.is.sexchrom)) {
            if (any(is_X)) {
                X_is_seXual <- any(is_Y) ||
                    any(is_nb_with_suffix) ||
                    any(is_nbxxx)
            } else {
                X_is_seXual <- TRUE  # or FALSE, won't make any difference
            }
        } else {
            X_is_seXual <- X.is.sexchrom
        }
        is_seXual <- is_X & X_is_seXual
        is_roman <- isRoman(sgsuffix) & !is_seXual
        is_U <- sgsuffix == "U"
        is_MT <- sgsuffix == "MT"
        is_M <- sgsuffix == "M"
        is_Wxxx <- .hasPrefix(sgsuffix, "W") & !is_W & !is_roman
        is_Zxxx <- .hasPrefix(sgsuffix, "Z") & !is_Z & !is_roman
        is_Xxxx <- .hasPrefix(sgsuffix, "X") & !is_X & !is_roman
        is_Yxxx <- .hasPrefix(sgsuffix, "Y") & !is_Y & !is_roman
        is_Uxxx <- .hasPrefix(sgsuffix, "U") & !is_U & !is_roman
        is_MTxxx <- .hasPrefix(sgsuffix, "MT") & !is_MT & !is_roman
        is_Mxxx <- .hasPrefix(sgsuffix, "M") & !is_M &
            !is_MT & !is_MTxxx & !is_roman
        ## The groups below must define a partitioning of the current super
        ## group i.e. any element in 'sgsuffix' must fall in exactly 1 group.
        is_xxx <- !is_roman & !is_nb_with_suffix &
            !is_W & !is_Z & !is_seXual & !is_Y & !is_U & !is_M & !is_MT &
            !is_nbxxx &
            !is_Wxxx & !is_Zxxx  & !is_Xxxx & !is_Yxxx &
                                   !is_Uxxx & !is_Mxxx & !is_MTxxx
        ## Group (a).
        if (any(is_roman)) {
            gsuffix <- sgsuffix[is_roman]
            ints <- as.integer(utils:::.roman2numeric(gsuffix))
            makeAndAssignProvIds(sgidx[is_roman], ints=ints)
        }
        ## Group (b).
        if (any(is_nb_with_suffix)) {
            gsuffix <- sgsuffix[is_nb_with_suffix]
            isnb_idx <- which(is_nb[is_nb_with_suffix])
            isnbA_idx <- which(is_nbA[is_nb_with_suffix])
            isnba_idx <- which(is_nba[is_nb_with_suffix])
            isnbB_idx <- which(is_nbB[is_nb_with_suffix])
            isnbb_idx <- which(is_nbb[is_nb_with_suffix])
            isnbL_idx <- which(is_nbL[is_nb_with_suffix])
            isnbR_idx <- which(is_nbR[is_nb_with_suffix])
            nb_ints <- as.integer(gsuffix[isnb_idx])
            gsuffixA <- gsuffix[isnbA_idx]
            nbA_ints <- as.integer(substr(gsuffixA,
                                          start=1L,
                                          stop=nchar(gsuffixA)-1L))
            gsuffixa <- gsuffix[isnba_idx]
            nba_ints <- as.integer(substr(gsuffixa,
                                          start=1L,
                                          stop=nchar(gsuffixa)-1L))
            gsuffixB <- gsuffix[isnbB_idx]
            nbB_ints <- as.integer(substr(gsuffixB,
                                          start=1L,
                                          stop=nchar(gsuffixB)-1L))
            gsuffixb <- gsuffix[isnbb_idx]
            nbb_ints <- as.integer(substr(gsuffixb,
                                          start=1L,
                                          stop=nchar(gsuffixb)-1L))
            gsuffixL <- gsuffix[isnbL_idx]
            nbL_ints <- as.integer(substr(gsuffixL,
                                          start=1L,
                                          stop=nchar(gsuffixL)-1L))
            gsuffixR <- gsuffix[isnbR_idx]
            nbR_ints <- as.integer(substr(gsuffixR,
                                          start=1L,
                                          stop=nchar(gsuffixR)-1L))
            ints <- integer(length(gsuffix))
            ints[isnb_idx] <- 7L * nb_ints
            ints[isnbA_idx] <- 7L * nbA_ints + 1L
            ints[isnba_idx] <- 7L * nba_ints + 2L
            ints[isnbB_idx] <- 7L * nbB_ints + 3L
            ints[isnbb_idx] <- 7L * nbb_ints + 4L
            ints[isnbL_idx] <- 7L * nbL_ints + 5L
            ints[isnbR_idx] <- 7L * nbR_ints + 6L
            makeAndAssignProvIds(sgidx[is_nb_with_suffix], ints=ints)
        }
        ## Group (c).
        if (any(is_W))
            makeAndAssignProvIds(sgidx[is_W])
        ## Group (d).
        if (any(is_Z))
            makeAndAssignProvIds(sgidx[is_Z])
        ## Group (e).
        if (any(is_seXual))
            makeAndAssignProvIds(sgidx[is_seXual])
        ## Group (f).
        if (any(is_Y))
            makeAndAssignProvIds(sgidx[is_Y])
        ## Group (g).
        if (any(is_U))
            makeAndAssignProvIds(sgidx[is_U])
        ## Group (h).
        if (any(is_M))
            makeAndAssignProvIds(sgidx[is_M])
        ## Group (i).
        if (any(is_MT))
            makeAndAssignProvIds(sgidx[is_MT])
        ## Group (j).
        if (any(is_nbxxx)) {
            gsuffix <- sgsuffix[is_nbxxx]
            ints1 <- as.integer(.getNbPart(gsuffix))
            ints2 <- as.integer(factor(.getPostNbPart(gsuffix)))
            ints <- (max(ints2) + 1L) * ints1 + ints2
            makeAndAssignProvIds(sgidx[is_nbxxx], ints=ints)
        }
        ## Group (k).
        if (any(is_Wxxx)) {
            gsuffix <- sgsuffix[is_Wxxx]
            ints <- as.integer(factor(.dropPrefix(gsuffix, 1L)))
            makeAndAssignProvIds(sgidx[is_Wxxx], ints=ints)
        }
        ## Group (l).
        if (any(is_Zxxx)) {
            gsuffix <- sgsuffix[is_Zxxx]
            ints <- as.integer(factor(.dropPrefix(gsuffix, 1L)))
            makeAndAssignProvIds(sgidx[is_Zxxx], ints=ints)
        }
        ## Group (m).
        if (any(is_Xxxx)) {
            gsuffix <- sgsuffix[is_Xxxx]
            ints <- as.integer(factor(.dropPrefix(gsuffix, 1L)))
            makeAndAssignProvIds(sgidx[is_Xxxx], ints=ints)
        }
        ## Group (n).
        if (any(is_Yxxx)) {
            gsuffix <- sgsuffix[is_Yxxx]
            ints <- as.integer(factor(.dropPrefix(gsuffix, 1L)))
            makeAndAssignProvIds(sgidx[is_Yxxx], ints=ints)
        }
        ## Group (o).
        if (any(is_Uxxx)) {
            gsuffix <- sgsuffix[is_Uxxx]
            ints <- as.integer(factor(.dropPrefix(gsuffix, 1L)))
            makeAndAssignProvIds(sgidx[is_Uxxx], ints=ints)
        }
        ## Group (p).
        if (any(is_Mxxx)) {
            gsuffix <- sgsuffix[is_Mxxx]
            ints <- as.integer(factor(.dropPrefix(gsuffix, 1L)))
            makeAndAssignProvIds(sgidx[is_Mxxx], ints=ints)
        }
        ## Group (q).
        if (any(is_MTxxx)) {
            gsuffix <- sgsuffix[is_MTxxx]
            ints <- as.integer(factor(.dropPrefix(gsuffix, 2L)))
            makeAndAssignProvIds(sgidx[is_MTxxx], ints=ints)
        }
        ## Group (r).
        if (any(is_xxx)) {
            gsuffix <- sgsuffix[is_xxx]
            ints <- as.integer(factor(gsuffix))
            makeAndAssignProvIds(sgidx[is_xxx], ints=ints)
        }
    }
    
    ## Longest prefixes first.
    assignProvIdsForSuperGroup(seqlevels, "CHR")
    assignProvIdsForSuperGroup(seqlevels, "chr")
    assignProvIdsForSuperGroup(seqlevels, "CH")
    assignProvIdsForSuperGroup(seqlevels, "ch")
    assignProvIdsForSuperGroup(seqlevels, "")
    
    ## Seqlevel ids.
    seqlevel_ids <- integer(length(prov_ids))
    oo <- order(prov_ids)
    seqlevel_ids[oo] <- seq_len(length(seqlevel_ids))
    
    ## Seqname ids.
    seqlevel_ids[as.integer(seqnames)]
}

orderSeqlevels <- 
  function(seqnames, X.is.sexchrom=NA)
  {
    if (missing(seqnames))
      seqnames <- character()
    order(rankSeqlevels(seqnames, X.is.sexchrom))
  }

