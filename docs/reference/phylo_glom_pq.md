# Agglomerate taxa closer than a cophenetic distance threshold

Merge taxa (ASV, OTU, ...) whose patristic (cophenetic) distance on the
`phy_tree` of a
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object is below a threshold `h`. Taxa are clustered with
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) on the
cophenetic distance matrix, cut at height `h` with
[`stats::cutree()`](https://rdrr.io/r/stats/cutree.html), and merged
with
[`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html).

This is a faster and more transparent alternative to
[`phyloseq::tip_glom()`](https://rdrr.io/pkg/phyloseq/man/tip_glom.html):
the clustering method is explicit, the mapping from taxa to clusters can
be returned instead of the merged object, and the number of merged taxa
is reported.

## Usage

``` r
phylo_glom_pq(
  physeq,
  h,
  hclust_method = "average",
  tax_adjust = 1L,
  rank_propagation = TRUE,
  return_map = FALSE,
  verbose = FALSE
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object with a non-empty `phy_tree` slot carrying branch lengths.

- h:

  (required) Numeric, the cophenetic distance threshold at which the
  tree-based dendrogram is cut. Taxa within `h` of each other are
  merged.

- hclust_method:

  Agglomeration method passed on to
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html). Default to
  `"average"` (UPGMA), as in
  [`phyloseq::tip_glom()`](https://rdrr.io/pkg/phyloseq/man/tip_glom.html).

- tax_adjust:

  Handling of taxonomic disagreements within a cluster. See
  [`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html).
  Default to 1 (phyloseq-compatible).

- rank_propagation:

  Logical, default TRUE, whether bad ranks are propagated to lower
  ranks. See
  [`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html).

- return_map:

  Logical, if TRUE return a data.frame mapping each taxon to its cluster
  instead of the agglomerated phyloseq object. Default to FALSE.

- verbose:

  Logical, if TRUE report the number of taxa before and after
  agglomeration. Default to FALSE.

## Value

A
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with agglomerated taxa, or a data.frame with columns `taxa` and
`cluster` when `return_map = TRUE`.

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## See also

[`phylo_glom_scan_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_scan_pq.md),
[`add_tree_pq()`](https://adrientaudiere.github.io/phylopq/reference/add_tree_pq.md),
[`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html),
[`phyloseq::tip_glom()`](https://rdrr.io/pkg/phyloseq/man/tip_glom.html)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
library(MiscMetabar)
data(data_fungi_mini)

pq <- add_tree_pq(data_fungi_mini, compute_brlen = TRUE)
phyloseq::ntaxa(pq)
#> [1] 45

pq_glom <- phylo_glom_pq(pq, h = 0.2, verbose = TRUE)
#> ✔ 45 taxa agglomerated into 26 taxa at h = 0.2.
phyloseq::ntaxa(pq_glom)
#> [1] 26

# Inspect the taxa-to-cluster mapping without merging
map <- phylo_glom_pq(pq, h = 0.2, return_map = TRUE)
head(map)
#>    taxa cluster
#> 1  ASV7       1
#> 2  ASV8       2
#> 3 ASV18       2
#> 4 ASV94       2
#> 5 ASV26       2
#> 6 ASV93       3
# }

if (FALSE) { # \dontrun{
# Keep the taxonomy strictly conservative across merged taxa
pq_glom <- phylo_glom_pq(pq, h = 0.2, tax_adjust = 2)
} # }
```
