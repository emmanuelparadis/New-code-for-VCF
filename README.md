# New code for reading VCF files

**UPDATE**: all code has now been included in [pegas](https://github.com/emmanuelparadis/pegas).

See below for some examples.

These functions provide new tools to handle VCF files. The new user-level functions are:

- `VCFheader` extracts the header of a VCF file as a single character string (can be printed in a more friendly way with `cat`).

- `VCFlabels` extracts the labels of the individuals.

- `VCFloci` extracts the information of the loci, by default all nine fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) are read. This returns an object of class "VCFinfo" with a `print` method.

- `is.snp` tests whether a locus is a SNP (this is now a generic function).

- `rangePOS` selects loci within a given range of POSition using an object of class "VCFinfo".

- `selectQUAL` selects loci above a given QUALity using an object of class "VCFinfo".

- `getINFO` extracts specific information from the INFO field (by default "DP").

The function `read.vcf` (already in pegas) has been rewritten with a new interface:

```r
read.vcf(file, from = 1, to = 1e4, which.loci = NULL, quiet = FALSE)
```

By default the first 10,000 loci. An alternative is to use `which.loci` which specifies which loci to read, typically as an output from the above function.

Another new feature is that `read.vcf` (and `VCFloci` too) can read both compressed (*.gz) or uncompressed files.

#### Example (with the chromosome Y from 1000 Genomes):

```r
> info <- VCFlociinfo("chrY.vcf")
Scanning file chrY.vcf 
 171.6615 / 171.6615 Mb
Done.
```

The INFO field is big so we remove it to get a nicer print:

```r
> INFO <- info$INFO
> info$INFO <- NULL
> info
      CHROM      POS         ID REF ALT QUAL FILTER FORMAT
1         Y  2655180 rs11575897   G   A  100   PASS     GT
2         Y  2655471          .   A   C  100   PASS     GT
3         Y  2655754          .   A   T  100   PASS     GT
4         Y  2655800          .   A   G  100   PASS     GT
....
.....
62039     Y 28770656          .   A   G  100   PASS     GT
62040     Y 28770756          .   C   G  100   PASS     GT
62041     Y 28770875          .   C   G  100   PASS     GT
62042     Y 28770931          .   T   C  100   PASS     GT
```

We can now identify the loci which are real SNPs (i.e., substitutions with only 2 alleles):

```r
> SNP <- is.snp(info)
> table(SNP)
SNP
FALSE  TRUE 
 1537 60505
```

We now read the SNPs:

```r
> X <- read.vcf("chrY.vcf", which.loci = which(SNP))
Reading 60505 loci.
Done
> X
Allelic data frame: 1233 individuals
                    60505 loci
```

And it's a simple matter to read the other (non-SNP) loci:

```r
> Y <- read.vcf("chrY.vcf", which.loci = which(!SNP))
Reading 1537 loci.
Done
> Y
Allelic data frame: 1233 individuals
                    1537 loci
```

A lot of things can be done from this (as usual with R). For instance, the help page `?read.vcf` gives code to draw the distribution of mutations along chromosome Y for non-SNP mutations (in red) and SNPs (in blue) for the whole chromosome (first plot) and on a restricted portion (second plot) marked with a dashed square:

![alt text](https://github.com/emmanuelparadis/New-code-for-VCF/blob/master/chrY.png)

