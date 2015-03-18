#include <R.h>
#include <Rinternals.h>

/* reads a binary stream of 'SIZE' bytes from file 'FILENAME'
   after skipping 'SKIP' bytes */
SEXP read_bin_pegas(SEXP FILENAME, SEXP SIZE, SEXP SKIP)
{
    SEXP res;
    const char *filename;
    FILE *fl;
    int n, i, sz;
    double skip;
    unsigned char *p;

    PROTECT(FILENAME = coerceVector(FILENAME, STRSXP));
    PROTECT(SIZE = coerceVector(SIZE, INTSXP)); // OK to have INT cause biggest chunk is 1e9
    PROTECT(SKIP = coerceVector(SKIP, REALSXP)); // !!! MUST be REAL cause file size can be > 2 Gb
    filename = CHAR(STRING_ELT(FILENAME, 0));
    sz = INTEGER(SIZE)[0];
    skip = REAL(SKIP)[0];
    PROTECT(res = allocVector(RAWSXP, sz));
    p = RAW(res);
    fl = fopen(filename, "r");
    fseek(fl, (long)skip, SEEK_SET);
    n = fread(p, 1, sz, fl);
    fclose(fl);
    UNPROTECT(4);
    return res;
}

/* change bytes into an integer */
static int raw2int(unsigned char *x, int a, int b)
{
	int i, k = 1, ans = 0;

	for (i = b; i >= a; i--, k *= 10)
		ans += ((int)x[i] - 48) * k;

	return ans;
}

void extract_substring(unsigned char *x, int a, int b, char *y)
{
	int i, j;

	for (i = a, j = 0; i <= b; i++, j++) y[j] = x[i];

	y[j] = '\0';
}

/* find the locations of end-of-lines (= linefeeds; LF) in a strem of bytes 'x'
   after skipping 'SKIP' bytes; once an LF is found, 'HOP' bytes are skipped */
SEXP findEOL_C(SEXP x, SEXP SKIP, SEXP HOP)
{
    int n, i, j, *p, *buf, nEOL, hop;
    unsigned char *xr;
    SEXP res;

    PROTECT(x = coerceVector(x, RAWSXP));
    PROTECT(SKIP = coerceVector(SKIP, INTSXP));
    PROTECT(HOP = coerceVector(HOP, INTSXP));
    n = LENGTH(x);
    xr = RAW(x);
    hop = INTEGER(HOP)[0];

    buf = (int*)R_alloc(n/hop, sizeof(int));

    i = INTEGER(SKIP)[0];
    nEOL = j = 0;
    while (i < n) {
	if (xr[i] == 0x0a) {
	    nEOL++;
	    buf[j] = i + 1;
	    j++;
	    i += hop;
	}
	i++;
    }

    PROTECT(res = allocVector(INTSXP, nEOL));
    p = INTEGER(res);

    for (i = 0; i < nEOL; i++) p[i] = buf[i];

    UNPROTECT(4);
    return res;
}

/* extract the (int) field of a VCF file given a stream
   of bytes 'x' with end-of-lines stored in 'EOL' after
   skipping 'nTABtoSKIP' TABs from the start of each line
   nTABtoSKIP = 1 -> POS
   nTABtoSKIP = 5 -> QUAL */
SEXP extract_POS(SEXP x, SEXP EOL, SEXP nTABtoSKIP)
{
    int n, i, j, k, a, *p, *eol;
    unsigned char *xr;
    SEXP res;

    PROTECT(x = coerceVector(x, RAWSXP));
    PROTECT(EOL = coerceVector(EOL, INTSXP));
    PROTECT(nTABtoSKIP = coerceVector(nTABtoSKIP, INTSXP));
    xr = RAW(x);
    n = LENGTH(EOL) - 1;
    eol = INTEGER(EOL);

    PROTECT(res = allocVector(INTSXP, n));
    p = INTEGER(res);

    for (i = 0; i < n; i++) {
	j = eol[i];
	for (k = 1; k <= INTEGER(nTABtoSKIP)[0]; k++) {
	    while (xr[j] != 0x09) j++;
	    j++;
	}
	a = j;
	while (xr[j] != 0x09) j++;
	p[i] = raw2int(xr, a, j - 1);
    }

    UNPROTECT(4);
    return res;
}

/* extract the (char) field of a VCF file given a stream
   of bytes 'x' with end-of-lines stored in 'EOL' after
   skipping 'nTABtoSKIP' TABs from the start of each line
   nTABtoSKIP = 0 -> CHROM
   nTABtoSKIP = 2 -> ID
   nTABtoSKIP = 3 -> REF
   nTABtoSKIP = 4 -> ALT
   nTABtoSKIP = 6 -> FILTER
   nTABtoSKIP = 7 -> INFO
   nTABtoSKIP = 8 -> FORMAT */
SEXP extract_REF(SEXP x, SEXP EOL, SEXP nTABtoSKIP)
{
    int n, i, j, k, a, *eol;
    unsigned char *xr;
    char str[1000];
    SEXP res;

    PROTECT(x = coerceVector(x, RAWSXP));
    PROTECT(EOL = coerceVector(EOL, INTSXP));
    PROTECT(nTABtoSKIP = coerceVector(nTABtoSKIP, INTSXP));
    xr = RAW(x);
    n = LENGTH(EOL) - 1;
    eol = INTEGER(EOL);

    PROTECT(res = allocVector(STRSXP, n));

    for (i = 0; i < n; i++) {
	j = eol[i];
	for (k = 1; k <= INTEGER(nTABtoSKIP)[0]; k++) {
	    while (xr[j] != 0x09) j++;
	    j++;
	}
	a = j;
	while (xr[j] != 0x09) j++;
	extract_substring(xr, a, j - 1, str);
	SET_STRING_ELT(res, i, mkChar(str));
    }

    UNPROTECT(4);
    return res;
}
