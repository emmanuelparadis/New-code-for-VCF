library(pegas)
source("getVCFinfo.R")

read.vcf <- function(file, from = 1, to = 1e4, which.loci = NULL, quiet = FALSE)
{
    f <- .VCFconnection(file)
    GZ <- if (inherits(f, "connection")) TRUE else FALSE

    if (is.null(which.loci)) which.loci <- from:to
    nLoci <- length(which.loci)

    meta <- .getMETAvcf(readBin(f, "raw", 1e5))
    labs <- strsplit(meta$LABELS, "\t")[[1]]
    nCol <- length(labs)
    n <- nCol - 9L
    hop <- 2L * nCol - 1L

    cache <- ls(env = .cacheVCF, all.names = TRUE)
    if (! file %in% cache) {
        if (!quiet) cat("File apparently not yet accessed:\n")
        info <- VCFlociinfo(file, what = "POS", quiet = quiet)
    }
    cache <- get(file, env = .cacheVCF)

    nChunks <- nrow(cache)
    obj <- vector("list", nLoci)
    locnms <- character(nLoci)

    if (GZ) open(f)

    ii <- 0L # number of loci read
    for (k in seq_len(nChunks)) {
        sel <- match(which.loci, cache$FROM[k]:cache$TO[k])
        sel <- sel[!is.na(sel)]
        if (!length(sel)) next

        if (k == 1) {
            ck <- cache$CHUNCK.SIZES[1L]
            skip <- 0L
        } else {
            ck <- cache$CHUNCK.SIZES[k]
            skip <- sum(cache$CHUNCK.SIZES[1L:(k - 1L)])
        }

        Y <- if (GZ) readBin(f, "raw", ck) else .Call("read_bin_pegas", file, ck, skip)
        skip <- if (k == 1) meta$position else 0L
        EOL <- .Call("findEOL_C", Y, skip, hop) # can multiply 'hop' by 2 if diploid

        for (i in sel) {
            start <- if (i == 1) skip + 1L else EOL[i - 1] + 1L
            end <- EOL[i] - 1L
            out <- .Call("build_factor_loci", Y[start:end], n)
            ii <- ii + 1L
            locnms[ii] <- out[[1L]]
            REF <- out[[2L]]
            ALT <- out[[3L]]
            geno <- out[[4L]]
            lv <- out[[5L]]

            ## substitute the allele names:
            tmp <- strsplit(ALT, ",")[[1L]]
            ## we start from last allele in case there are more than 10 alleles...
            for (j in length(tmp):1)
                lv <- gsub(as.character(j), tmp[j], lv)
            ## ... and we finish with the reference allele:
            lv <- gsub("0", REF, lv)

            attr(geno, "levels") <- lv
            class(geno) <- "factor"
            obj[[ii]] <- geno
            #if (!quiet && !(ii %% 100))
                cat("\rReading", ii, "/", nLoci, "loci")
            ##if (ii == 47193) browser()
        }
    }
    if (GZ) close(f)
    if (!quiet) cat("\rReading", ii, "/", nLoci, "loci.\nDone.\n")

    names(obj) <- locnms
    class(obj) <- c("loci", "data.frame")
    attr(obj, "locicol") <- 1:nLoci
    rownames(obj) <- labs[-(1:9)]
    obj
}
