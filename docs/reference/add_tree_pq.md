# Attach, repair and prune a phylogenetic tree onto a phyloseq object

Safely add a phylogenetic tree to the `phy_tree` slot of a
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object. The tree may be a `phylo` object, a path to a Newick or Nexus
file, or `NULL`, in which case a taxonomy-based tree is built on the fly
with
[`taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.md).

Tip labels and taxa names rarely match exactly. `add_tree_pq()`
reconciles them explicitly: tips absent from the phyloseq object are
dropped, taxa absent from the tree are pruned, and the function reports
what it removed instead of failing silently or letting `phyloseq`
intersect quietly.

## Usage

``` r
add_tree_pq(
  physeq,
  tree = NULL,
  drop_tips = TRUE,
  prune_taxa = TRUE,
  root = c("none", "midpoint", "outgroup"),
  outgroup = NULL,
  ladderize = TRUE,
  compute_brlen = FALSE,
  force = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object

- tree:

  The tree to attach. Either a `phylo` object, a length-one character
  giving the path to a Newick (`.nwk`, `.tre`, `.newick`) or Nexus
  (`.nex`, `.nexus`) file, or `NULL` (default) to build a taxonomic tree
  with
  [`taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.md).

- drop_tips:

  Logical, if TRUE (default) tips of `tree` that are not taxa of
  `physeq` are dropped with
  [`ape::drop.tip()`](https://rdrr.io/pkg/ape/man/drop.tip.html).

- prune_taxa:

  Logical, if TRUE (default) taxa of `physeq` that are not tips of
  `tree` are pruned with
  [`phyloseq::prune_taxa()`](https://rdrr.io/pkg/phyloseq/man/prune_taxa-methods.html).

- root:

  One of `"none"` (default, keep the tree as it is), `"midpoint"`
  (midpoint rooting, requires the `phangorn` package) or `"outgroup"`
  (root on `outgroup` with
  [`ape::root()`](https://rdrr.io/pkg/ape/man/root.html)).

- outgroup:

  Character vector of tip labels used as outgroup. Only used when
  `root = "outgroup"`.

- ladderize:

  Logical, if TRUE (default) the tree is ladderized with
  [`ape::ladderize()`](https://rdrr.io/pkg/ape/man/ladderize.html).

- compute_brlen:

  Logical, if TRUE and the tree has no branch lengths, arbitrary branch
  lengths are computed with
  [`ape::compute.brlen()`](https://rdrr.io/pkg/ape/man/compute.brlen.html).
  This is required by downstream distance-based tools (e.g.
  [`phylo_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_pq.md)
  or
  [`phyloseq::UniFrac()`](https://rdrr.io/pkg/phyloseq/man/UniFrac-methods.html))
  but the resulting lengths carry no evolutionary meaning. Default to
  FALSE.

- force:

  Logical, if TRUE an existing `phy_tree` slot is replaced. Default to
  FALSE, in which case a non-empty slot triggers an error.

- verbose:

  Logical, if TRUE report the number of dropped tips and pruned taxa.
  Default to FALSE.

- ...:

  Additional arguments passed on to
  [`taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.md)
  when `tree` is NULL.

## Value

A
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with a `phy_tree` slot whose tip labels match its taxa names
exactly.

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## See also

[`taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.md),
[`phylo_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_pq.md),
[`phyloseq::phy_tree()`](https://rdrr.io/pkg/phyloseq/man/phy_tree-methods.html)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
library(MiscMetabar)
#> Loading required package: ggplot2
#> Loading required package: dplyr
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
data(data_fungi_mini)

# Build a taxonomic tree and attach it in one call
pq_tree <- add_tree_pq(data_fungi_mini, verbose = TRUE)
#> ✔ Attached a tree with 45 tips.
#> ℹ 0 tips dropped from the tree.
#> ℹ 0 taxa pruned from the phyloseq object.
pq_tree
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 45 taxa and 137 samples ]
#> sample_data() Sample Data:       [ 137 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 45 taxa by 12 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 45 tips and 80 internal nodes ]
#> refseq()      DNAStringSet:      [ 45 reference sequences ]

# Attach an existing phylo object, adding arbitrary branch lengths
tree <- taxo2tree(data_fungi_mini)
pq_brlen <- add_tree_pq(data_fungi_mini, tree, compute_brlen = TRUE)
ape::is.ultrametric(phyloseq::phy_tree(pq_brlen))
#> [1] TRUE

# A tree covering only part of the taxa prunes the phyloseq object
small_tree <- ape::drop.tip(tree, tree$tip.label[1:10])
pq_small <- add_tree_pq(data_fungi_mini, small_tree, verbose = TRUE)
#> ✔ Attached a tree with 35 tips.
#> ℹ 0 tips dropped from the tree.
#> ℹ 10 taxa pruned from the phyloseq object.
phyloseq::ntaxa(pq_small)
#> [1] 35
# }

if (FALSE) { # \dontrun{
# Read a tree from a Newick file written by an external program
pq <- add_tree_pq(data_fungi_mini, "my_tree.nwk", root = "midpoint")

# Replace an existing tree
pq <- add_tree_pq(pq, tree, force = TRUE)
} # }
```
