# Scan several cophenetic thresholds before agglomerating

Compute the number of taxa that
[`phylo_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_pq.md)
would return for a range of cophenetic thresholds, without building the
agglomerated objects. Useful to choose an informed value of `h` before
running
[`phylo_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_pq.md).

## Usage

``` r
phylo_glom_scan_pq(physeq, h_values, hclust_method = "average")
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object with a non-empty `phy_tree` slot carrying branch lengths.

- h_values:

  (required) Numeric vector of cophenetic distance thresholds to scan.

- hclust_method:

  Agglomeration method passed on to
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html). Default to
  `"average"` (UPGMA), as in
  [`phyloseq::tip_glom()`](https://rdrr.io/pkg/phyloseq/man/tip_glom.html).

## Value

A data.frame with one row per value of `h_values` and the columns `h`
(the threshold), `n_taxa` (the resulting number of taxa) and `prop_taxa`
(that number divided by the initial number of taxa).

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## See also

[`phylo_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_pq.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
library(MiscMetabar)
data(data_fungi_mini)

pq <- add_tree_pq(data_fungi_mini, compute_brlen = TRUE)
phylo_glom_scan_pq(pq, h_values = seq(0, 1, by = 0.1))
#>      h n_taxa prop_taxa
#> 1  0.0     45 1.0000000
#> 2  0.1     36 0.8000000
#> 3  0.2     26 0.5777778
#> 4  0.3     19 0.4222222
#> 5  0.4     16 0.3555556
#> 6  0.5     11 0.2444444
#> 7  0.6     11 0.2444444
#> 8  0.7     11 0.2444444
#> 9  0.8     11 0.2444444
#> 10 0.9     11 0.2444444
#> 11 1.0     11 0.2444444
# }

if (FALSE) { # \dontrun{
scan <- phylo_glom_scan_pq(pq, h_values = seq(0, 1, by = 0.01))
plot(scan$h, scan$n_taxa, type = "l")
} # }
```
