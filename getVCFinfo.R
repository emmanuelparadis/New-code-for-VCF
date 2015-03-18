dyn.load("readVCFbin.so")

.VCFconnection <- function(file)
{
    if (length(grep("\\.vcf$", file))) return(file)
    if (length(grep("\\.vcf.gz$", file))) return(gzcon(gzfile(file)))
    stop("file name does not end with '.vcf' or '.vcf.gz'")
}

.getMETAvcf <- function(x, position.only = FALSE)
{
    i <- 1L
    while (!identical(x[i + 0:5], charToRaw("#CHROM"))) i <- i + 1L
    j <- i + 6L
    while (x[j] != charToRaw("\n")) j <- j + 1L
    if (position.only) return(j)
    list(HEADER = rawToChar(x[1:(i - 1)]),
         LABELS = rawToChar(x[i:(j - 1)]), # avoid returning the last "\n"
         position = j)
}

VCFheader <- function(file)
{
    f <- .VCFconnection(file)
    x <- readBin(f, "raw", 1e5)
    .getMETAvcf(x)$HEADER
}

VCFlabels <- function(file)
{
    f <- .VCFconnection(file)
    x <- readBin(f, "raw", 1e5)
    strsplit(.getMETAvcf(x)$LABELS, "\t")[[1]]
}

chr1 <- "../Téléchargements/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
chr22 <- "../Téléchargements/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
titi <- "/media/paradis/12674652-3312-4cf3-9e96-3b970feca8c6/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"
fl <- "tmp.vcf.gz"
fl2 <- "tmp.vcf"

VCFlociinfo <- function(file, what = "POS", chunck.size = 1e9)
{
    f <- .VCFconnection(file)
    GZ <- if (inherits(f, "connection")) TRUE else FALSE

    FIELDS <- c("CHROM", "POS", "ID", "REF", "ALT",
                "QUAL", "FILTER", "INFO", "FORMAT")

    what <-
        if (identical(what, "all")) 1:9 else match(what, FIELDS)

    obj <- vector("list", 9L)

    if (GZ) {
        open(f)
    } else {
        sz <- file.info(file)$size
        if (is.na(sz))
            stop(paste("cannot find information on file", sQuote(file)))
        left.to.scan <- sz
        scanned <- 0 # !!! NOT integer please ;)
    }

    headerFound <- FALSE
    cat("Scanning file", file, "\n")
    ncycle <- 0L

    repeat {
        if (!GZ) {
            if (left.to.scan > chunck.size) {
                ck <- chunck.size
                left.to.scan <- left.to.scan - ck
            } else {
                ck <- left.to.scan
                left.to.scan <- 0L
            }
            Y <- .Call("read_bin_pegas", file, ck, scanned)
            scanned <- scanned + ck
        } else {
            Y <- readBin(f, "raw", chunck.size)
            if (!length(Y)) break
        }

        if (GZ) cat("\r", ncycle*1e3, "Mb")
        else cat("\r", ncycle*1e3, "/", sz/1e6, "Mb")

        ncycle <- ncycle + 1L

        if (!headerFound) {
            meta <- .getMETAvcf(Y)
            skip <- meta$position - 3L
            nCol <- length(gregexpr("\t", meta$LABELS)[[1]]) + 1L
            headerFound <- TRUE
        } else skip <- 0L

        hop <- 2L * nCol - 1L
        EOL <- .Call("findEOL_C", Y, skip, hop) # can multiply 'hop' by 2 if diploid

        ## we assume there are no extra blank line!

        for (i in what) {
            lib <- if (i %in% c(2, 6)) "extract_POS" else "extract_REF"
            tmp <- .Call(lib, Y, EOL, i - 1L)
            obj[[i]] <- c(obj[[i]], tmp)
        }

        if (!GZ && !left.to.scan) break
    }

    if (GZ) close(f)
    else cat("\r", scanned/1e6, "/", sz/1e6, "Mb")
    cat("\nDone.\n")

    names(obj) <- FIELDS
    obj <- obj[!sapply(obj, is.null)]
    as.data.frame(obj, stringsAsFactors = FALSE)
}
