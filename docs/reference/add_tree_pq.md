# Attach, repair and prune a phylogenetic tree onto a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Safely add a phylogenetic tree to the `phy_tree` slot of a
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object. By default the tree is supplied by the user as a `phylo` object
(`ape` package) or as a path to a Newick or Nexus file. Set
`use_taxo_to_build_tree = TRUE` to build a taxonomy-based tree on the
fly with
[`taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.md)
instead.

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
  use_taxo_to_build_tree = FALSE,
  drop_tips = TRUE,
  prune_taxa = FALSE,
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

  (required unless `use_taxo_to_build_tree` is TRUE) The tree to attach.
  Either a `phylo` object (`ape` package) or a length-one character
  giving the path to a Newick (`.nwk`, `.tre`, `.newick`) or Nexus
  (`.nex`, `.nexus`) file.

- use_taxo_to_build_tree:

  Logical, if TRUE a taxonomy-based tree is built from the `tax_table`
  slot with
  [`taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.md)
  and `tree` must be left NULL. Default to FALSE.

- drop_tips:

  Logical, if TRUE (default) tips of `tree` that are not taxa of
  `physeq` are dropped with
  [`ape::drop.tip()`](https://rdrr.io/pkg/ape/man/drop.tip.html).

- prune_taxa:

  Logical, if TRUE taxa of `physeq` that are not tips of `tree` are
  pruned with
  [`phyloseq::prune_taxa()`](https://rdrr.io/pkg/phyloseq/man/prune_taxa-methods.html).
  Default to FALSE, in which case such taxa trigger an error: dropping
  data is opt-in, since a tree that does not cover every taxon usually
  means something went wrong upstream.

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
  when `use_taxo_to_build_tree` is TRUE.

## Value

A
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with a `phy_tree` slot whose tip labels match its taxa names
exactly.

## See also

[`taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.md),
[`phylo_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_pq.md),
[`MiscMetabar::build_phytree_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/build_phytree_pq.html),
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

# Attach a phylo object
tree <- taxo2tree(data_fungi_mini)
pq_tree <- add_tree_pq(data_fungi_mini, tree, verbose = TRUE)
#> ✔ Attached a tree with 45 tips.
#> ℹ 0 tips dropped from the tree.
#> ℹ 0 taxa pruned from the phyloseq object.
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
pq_tree
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 45 taxa and 137 samples ]
#> sample_data() Sample Data:       [ 137 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 45 taxa by 12 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 45 tips and 80 internal nodes ]
#> refseq()      DNAStringSet:      [ 45 reference sequences ]

# Build a taxonomic tree and attach it in one call, adding arbitrary
# branch lengths
pq_brlen <- add_tree_pq(
  data_fungi_mini,
  use_taxo_to_build_tree = TRUE,
  compute_brlen = TRUE
)
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
ape::is.ultrametric(phyloseq::phy_tree(pq_brlen))
#> [1] TRUE

# A tree covering only part of the taxa errors unless pruning is opt-in
small_tree <- ape::drop.tip(tree, tree$tip.label[1:10])
try(add_tree_pq(data_fungi_mini, small_tree))
#> Error in add_tree_pq(data_fungi_mini, small_tree) : 
#>   10 taxa of `physeq` are absent from `tree`.
#> ℹ Use `prune_taxa = TRUE` to prune them.
pq_small <- add_tree_pq(
  data_fungi_mini,
  small_tree,
  prune_taxa = TRUE,
  verbose = TRUE
)
#> ✔ Attached a tree with 35 tips.
#> ℹ 0 tips dropped from the tree.
#> ℹ 10 taxa pruned from the phyloseq object.
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
phyloseq::ntaxa(pq_small)
#> [1] 35
# }

if (FALSE) { # \dontrun{
# Read a tree from a Newick file written by an external program
pq <- add_tree_pq(data_fungi_mini, "my_tree.nwk", root = "midpoint")

# Replace an existing tree
pq <- add_tree_pq(pq, tree, force = TRUE)

# Attach a real phylogeny inferred from the refseq slot with
# MiscMetabar::build_phytree_pq(), which returns a list of trees. Such a
# tree already carries branch lengths, so `compute_brlen` is not needed and
# the result can feed distance-based tools directly.
set.seed(22)
df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 9000)
phytree <- MiscMetabar::build_phytree_pq(df)

pq_ml <- add_tree_pq(df, phytree$ML$tree, root = "midpoint")
phylo_glom_pq(pq_ml, h = 0.05, verbose = TRUE)

# The other trees of the list are attached the same way
pq_upgma <- add_tree_pq(df, phytree$UPGMA)
pq_nj <- add_tree_pq(df, phytree$NJ, root = "midpoint")
} # }
```
