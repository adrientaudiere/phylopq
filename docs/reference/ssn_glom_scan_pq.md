# Scan several similarity thresholds before agglomerating

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Compute the number of network sequence clusters that
[`ssn_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/ssn_glom_pq.md)
would return for a range of similarity thresholds, without building the
agglomerated objects. Useful to choose an informed value of `id`, as a
single-linkage network is sensitive to it.

vsearch is run **once**, at the lowest threshold of `id_values`, and the
resulting pairwise identities are then filtered at each threshold in
turn — so scanning twenty thresholds costs one alignment run, not
twenty.

## Usage

``` r
ssn_glom_scan_pq(
  physeq,
  id_values,
  iddef = 1,
  vsearchpath = NULL,
  vsearch_args = "",
  keep_temporary_files = FALSE,
  verbose = FALSE
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object with a non-empty `refseq` slot.

- id_values:

  (required) Numeric vector of similarity thresholds in \[0, 1\] to
  scan.

- iddef:

  Integer, the identity definition used by vsearch: 0 (CD-HIT
  definition), 1 (edit distance, the default here and the definition
  used by Forster et al. 2019), 2 (edit distance excluding terminal
  gaps), 3 (marine-biological definition) or 4 (edit distance excluding
  gaps of any kind).

- vsearchpath:

  Path to the vsearch executable. Default to NULL, in which case it is
  looked up with
  [`MiscMetabar::find_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/reference/find_vsearch.html).

- vsearch_args:

  A length-one character vector of further arguments passed on to
  vsearch, e.g. `"--threads 4"`. Note that `--strand` is not a valid
  option of `--allpairs_global`. Default to `""`.

- keep_temporary_files:

  Logical, if TRUE the FASTA and `.uc` files handed to and produced by
  vsearch are kept. Default to FALSE.

- verbose:

  Logical, if TRUE report the number of edges, of clusters and of taxa
  before and after agglomeration. Default to FALSE.

## Value

A data.frame with one row per value of `id_values` and the columns `id`
(the threshold), `n_edges` (the number of pairs at or above it),
`n_taxa` (the resulting number of clusters) and `prop_taxa` (that number
divided by the initial number of taxa).

## See also

[`ssn_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/ssn_glom_pq.md),
[`phylo_glom_scan_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_scan_pq.md)

## Author

Adrien Taudière

## Examples

``` r
library(MiscMetabar)
data(data_fungi_mini)
df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 5000)
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 45
#> Number of filtered-out ASV 0
#> Number of kept ASV 45
#> Number of kept samples 137

ssn_glom_scan_pq(df, id_values = seq(0.8, 1, by = 0.05))
#>     id n_edges n_taxa prop_taxa
#> 1 0.80      56     21 0.4666667
#> 2 0.85      45     24 0.5333333
#> 3 0.90      36     28 0.6222222
#> 4 0.95      31     29 0.6444444
#> 5 1.00       0     45 1.0000000

if (FALSE) { # \dontrun{
scan <- ssn_glom_scan_pq(df, id_values = seq(0.7, 1, by = 0.01))
plot(scan$id, scan$n_taxa, type = "l", xlab = "id", ylab = "Nb of NSC")
} # }
```
