# Species delimitation of a phyloseq object with ABGD or ASAP

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Run a single-locus species-delimitation analysis on the reference
sequences of a
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object and map the resulting partitions back onto its taxa.

Two barcode-gap methods are available through the `delimtools` package:
ABGD (Automatic Barcode Gap Discovery, Puillandre et al. 2012) and ASAP
(Assemble Species by Automatic Partitioning, Puillandre et al. 2021).
Both are external command-line programs: either give the path to the
executable with `exe`, or pass a result file downloaded from the
corresponding web server with `webserver`.

By default taxa belonging to the same partition are merged with
[`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html),
so `delim_pq()` acts as a post-clustering step comparable to
[`MiscMetabar::postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.html),
but driven by a barcode gap rather than by a fixed similarity threshold.

## Usage

``` r
delim_pq(
  physeq,
  method = c("asap", "abgd"),
  exe = NULL,
  model = 3,
  slope = 1.5,
  webserver = NULL,
  outfolder = NULL,
  align = TRUE,
  align_method = c("decipher", "mafft"),
  mafft_exec = NULL,
  merge_taxa = TRUE,
  tax_adjust = 1L,
  rank_propagation = TRUE,
  keep_temporary_files = FALSE,
  verbose = FALSE
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object with a non-empty `refseq` slot.

- method:

  One of `"asap"` (default) or `"abgd"`.

- exe:

  Path to the ABGD or ASAP executable. Ignored when `webserver` is
  given. Default to NULL, in which case the executable is looked up with
  [`is_delim_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_delim_installed.md):
  first the `phylopq.asappath` / `phylopq.abgdpath` option, then a copy
  installed by
  [`install_asap()`](https://adrientaudiere.github.io/phylopq/reference/install_asap.md)
  /
  [`install_abgd()`](https://adrientaudiere.github.io/phylopq/reference/install_abgd.md),
  then the system `PATH`.

- model:

  Integer, the evolutionary model used by the external program. 0:
  Kimura-2P, 1: Jukes-Cantor, 2: Tamura-Nei, 3: simple p-distance
  (default).

- slope:

  Numeric, relative gap width. Only used when `method = "abgd"`. Default
  to 1.5.

- webserver:

  A result file obtained from the ABGD (`.txt`) or ASAP (`.csv`) web
  server. When given, no executable is required. Default to NULL.

- outfolder:

  Path to the folder where the external program writes its output.
  Default to NULL, i.e. a temporary location.

- align:

  Logical, if TRUE (default) reference sequences of unequal length are
  aligned with
  [`align_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/align_pq.html)
  before being submitted. Both ABGD and ASAP require an aligned FASTA
  and silently return nothing otherwise, so set `align = FALSE` only
  when `refseq` is already aligned.

- align_method:

  Aligner used when `align = TRUE`, either `"decipher"` (default,
  pure R) or `"mafft"` (external program, much faster on large `refseq`
  slots). See
  [`align_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/align_pq.html).

- mafft_exec:

  Path to the MAFFT executable. Only used when `align_method = "mafft"`.
  Default to NULL, i.e. the usual lookup of
  [`is_mafft_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_mafft_installed.html).

- merge_taxa:

  Logical, if TRUE (default) taxa of the same partition are merged and a
  phyloseq object is returned. If FALSE, the partition table is returned
  untouched.

- tax_adjust:

  Handling of taxonomic disagreements within a partition. See
  [`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html).
  Default to 1.

- rank_propagation:

  Logical, default TRUE, whether bad ranks are propagated to lower
  ranks. See
  [`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html).

- keep_temporary_files:

  Logical, if TRUE the temporary FASTA file submitted to the external
  program is kept. Default to FALSE.

- verbose:

  Logical, if TRUE report the number of taxa before and after
  delimitation. Default to FALSE.

## Value

A
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object whose taxa are the delimited species, or a data.frame with
columns `taxa` and `partition` when `merge_taxa = FALSE`.

## Details

ABGD and ASAP are not available on Windows. Taxa that the external
program does not return a partition for are dropped by
[`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html)
and reported with a warning.

Both programs need an aligned FASTA. Metabarcoding reference sequences
are almost never aligned, hence the `align = TRUE` default; alignment of
a large `refseq` slot can be slow with the default `"decipher"` backend,
in which case `align_method = "mafft"` is usually much faster.

## References

Puillandre N., Lambert A., Brouillet S., Achaz G. (2012) ABGD, Automatic
Barcode Gap Discovery for primary species delimitation. *Molecular
Ecology* 21(8):1864-1877.
[doi:10.1111/j.1365-294X.2011.05239.x](https://doi.org/10.1111/j.1365-294X.2011.05239.x)

Puillandre N., Brouillet S., Achaz G. (2021) ASAP: assemble species by
automatic partitioning. *Molecular Ecology Resources* 21:609-620.
[doi:10.1111/1755-0998.13281](https://doi.org/10.1111/1755-0998.13281)

## See also

[`delim_multi_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_multi_pq.md),
[`is_delim_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_delim_installed.md),
[`align_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/align_pq.html),
[`MiscMetabar::postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.html),
[`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html)

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
data(data_fungi_mini)
pq_asap <- delim_pq(data_fungi_mini, method = "asap", verbose = TRUE, slope=0.5)
#> /*
#>  ASAP (Agglomerate Specimens by Automatic Processing)
#>  will delineate species in your dataset in a few moments.
#>  Remember that the final cut remains yours!
#> */
#> > asap is reading the fasta file and computing the distance matrix
#> done read
#> done mat
#> End of matrix distance
#>   45 input sequences
#> > asap is building and testing all partitions
#>   > asap has finished building and testing all partitions
#>   
#> > 10 Best asap scores (probabilities evaluated with seq length:423)
#>   distance  #species   #spec w/rec  p-value pente asap-score
#>    0.1274       27            27  1.557e-01 5.832780e-01     3.000000 
#> *  0.0268       31            31  2.375e-01 3.069919e-01     5.500000 
#>    0.1093       28            28  8.184e-01 1.126107e+00     6.000000 
#>    0.2025       21            21  4.504e-03 1.604537e-01     6.500000 
#>    0.1410       25            25  1.084e-02 2.166865e-01     6.500000 
#> *  0.0143       32            32  7.844e-01 3.382643e-01     6.500000 
#> *  0.0442       30            30  7.884e-01 3.041927e-01     8.000000 
#>    0.1482       24            24  9.501e-01 4.065041e-01     8.500000 
#> *  0.0086       33            33  7.784e-01 1.526592e-01     9.500000 
#>    0.1734       23            23  8.623e-01 2.993369e-01     10.000000 
#> > asap is creating text and graphical output
#> > results were write 
#>   partition results are logged in: 
#>  Full results: /tmp/RtmpSIJqWk/phylopq_delim_c350d37d4fc94/delim_pq_asap.f.all
#>  The rank below is given by asap_score
#>  The csv file of rank x: /tmp/RtmpSIJqWk/phylopq_delim_c350d37d4fc94/Partition_x
#>  The result file of rank x: /tmp/RtmpSIJqWk/phylopq_delim_c350d37d4fc94/Partition_x
#>  The graphic outputs are: /tmp/RtmpSIJqWk/phylopq_delim_c350d37d4fc94/*.svg
#>  XML spart file is: /tmp/RtmpSIJqWk/phylopq_delim_c350d37d4fc94/delim_pq_asap.f.spart.xml
#>  Spart file is: /tmp/RtmpSIJqWk/phylopq_delim_c350d37d4fc94/delim_pq_asap.f.spart
#> > asap computation times were:
#>    0m  0s to read file and compute distance
#>    0m  1s to compute and test all partitions
#>   --------------
#>    0m  1s total
#> ✔ 45 taxa delimited into 27 species by asap.
phyloseq::ntaxa(pq_asap)
#> [1] 27

if (FALSE) { # \dontrun{
# Give an explicit path to the executable and inspect the partitions only
partitions <- delim_pq(
  data_fungi_mini,
  method = "abgd",
  exe = "/usr/local/bin/abgd",
  slope = 0.5,
  merge_taxa = FALSE
)
head(partitions)

# Reuse a result file downloaded from the ASAP web server
pq <- delim_pq(data_fungi_mini, method = "asap", webserver = "asap.csv")

# Align with MAFFT instead of DECIPHER, much faster on large refseq slots
pq <- delim_pq(data_fungi_mini, method = "asap", align_method = "mafft")
} # }
```
