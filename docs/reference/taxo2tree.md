# Convert taxonomy dataframe to phylogenetic tree

Creates a phylo object from a taxonomy table with hierarchical taxonomic
ranks as columns.

## Usage

``` r
taxo2tree(
  physeq,
  ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  internal_node_singletons = TRUE,
  use_taxa_names = TRUE
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object

- ranks:

  Character vector specifying the column names to use as taxonomic
  ranks, ordered from highest to lowest. By default: c("Domain",
  "Phylum", "Class", "Order", "Family", "Genus", "Species").

- internal_node_singletons:

  Logical, if TRUE, create internal nodes for singleton. If FALSE,
  internal nodes with only one descendant are discarded.

- use_taxa_names:

  Logical, if TRUE (default), use the taxa names (rownames, e.g., ASV_1,
  ASV_2) as terminal leaves. If FALSE, collapse identical taxonomy paths
  and use the lowest rank value as tip labels. This is useful for
  creating cleaner trees that show only taxonomic structure without
  individual ASV/OTU names.

## Value

A phylo object (ape package) representing the taxonomic tree.

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

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
tree <- taxo2tree(data_fungi_mini,
  ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
)
tree
#> 
#> Phylogenetic tree with 45 tips and 56 internal nodes.
#> 
#> Tip labels:
#>   ASV7, ASV8, ASV18, ASV26, ASV94, ASV93, ...
#> Node labels:
#>   Fungi, Basidiomycota, Agaricomycetes, Russulales, Stereaceae, Stereum, ...
#> 
#> Rooted; no branch length.

# Without internal node singletons
tree_wo_singletons <- taxo2tree(data_fungi_mini,
  ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
  internal_node_singletons = FALSE
)

# Without taxa names (collapse identical paths)
tree_no_taxa <- taxo2tree(data_fungi_mini,
  ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
  use_taxa_names = FALSE
)
# }

if (FALSE) { # \dontrun{
# Attach the taxonomic tree to a phyloseq object and plot it with ggtree
library(MiscMetabar)
library(ggtree)
data(data_fungi_mini)
data_fungi_mini@phy_tree <- phyloseq::phy_tree(
  taxo2tree(data_fungi_mini,
    ranks = c(
      "Domain", "Phylum", "Class", "Order", "Family", "Genus",
      "Genus_species"
    )
  )
)

ggtree(data_fungi_mini@phy_tree) +
  geom_nodelab(size = 2, nudge_x = -0.2, nudge_y = 0.6) +
  geom_tiplab()
} # }
```
