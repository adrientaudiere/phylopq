# Species delimitation of a phyloseq object with ABGD or ASAP

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
  given. Default to NULL, in which case the executable is looked up on
  the `PATH` with
  [`is_delim_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_delim_installed.md).

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
  [`DECIPHER::AlignSeqs()`](https://rdrr.io/pkg/DECIPHER/man/AlignSeqs.html)
  before being submitted. Both ABGD and ASAP require an aligned FASTA
  and silently return nothing otherwise, so set `align = FALSE` only
  when `refseq` is already aligned.

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

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

ABGD and ASAP are not available on Windows. Taxa that the external
program does not return a partition for are dropped by
[`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html)
and reported with a warning.

Both programs need an aligned FASTA. Metabarcoding reference sequences
are almost never aligned, hence the `align = TRUE` default; alignment of
a large `refseq` slot can be slow.

## References

Puillandre N., Lambert A., Brouillet S., Achaz G. (2012) ABGD, Automatic
Barcode Gap Discovery for primary species delimitation. *Molecular
Ecology* 21(8):1864-1877.
[doi:10.1111/j.1365-294X.2011.05239.x](https://doi.org/10.1111/j.1365-294X.2011.05239.x)

Puillandre N., Brouillet S., Achaz G. (2021) ASAP: assemble species by
automatic partitioning. *Molecular Ecology Resources* 21:609-620.
[doi:10.1111/1755-0998.13281](https://doi.org/10.1111/1755-0998.13281)

## See also

[`is_delim_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_delim_installed.md),
[`MiscMetabar::postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.html),
[`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # phylopq::is_delim_installed("asap")
library(MiscMetabar)
data(data_fungi_mini)
pq_asap <- delim_pq(data_fungi_mini, method = "asap", verbose = TRUE)
phyloseq::ntaxa(pq_asap)

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
} # }
}
```
