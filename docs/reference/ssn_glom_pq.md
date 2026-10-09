# Agglomerate taxa into network sequence clusters (NSC)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Group the taxa of a
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object into *network sequence clusters* (NSC), following the
sequence-similarity network (SSN) approach of Forster et al. (2019). All
pairs of reference sequences are compared with
`vsearch --allpairs_global`, every pair at least `id` similar becomes an
edge of an undirected graph, and each connected component of that graph
is one cluster. Taxa of the same component are then merged with
[`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html).

A connected component is a *single-linkage* group: two sequences end up
in the same NSC as soon as a chain of pairwise similarities links them,
even when they are themselves below `id`. This is the point of the
method — it recovers the loose, chained variation of a species that a
fixed centroid threshold splits apart — but it also means an NSC can be
much wider than `id`, and that a single bridging sequence can fuse two
clusters. Use
[`ssn_glom_scan_pq()`](https://adrientaudiere.github.io/phylopq/reference/ssn_glom_scan_pq.md)
to see how the number of clusters responds to `id` before committing to
a value.

Unlike
[`MiscMetabar::postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.html),
which assigns each sequence to a centroid, and unlike
[`phylo_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_pq.md),
which cuts a tree, `ssn_glom_pq()` needs no centroid and no tree: only
the `refseq` slot.

## Usage

``` r
ssn_glom_pq(
  physeq,
  id = 0.97,
  iddef = 1,
  vsearchpath = NULL,
  vsearch_args = "",
  tax_adjust = 1L,
  rank_propagation = TRUE,
  return_map = FALSE,
  return_graph = FALSE,
  keep_temporary_files = FALSE,
  verbose = FALSE
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object with a non-empty `refseq` slot.

- id:

  Numeric in \[0, 1\], the pairwise-similarity threshold above which two
  sequences are linked by an edge. Default to 0.97.

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

- return_graph:

  Logical, if TRUE return the `igraph` sequence-similarity network
  itself, with the cluster membership stored as the `nsc` vertex
  attribute. Default to FALSE. Takes precedence over `return_map`.

- keep_temporary_files:

  Logical, if TRUE the FASTA and `.uc` files handed to and produced by
  vsearch are kept. Default to FALSE.

- verbose:

  Logical, if TRUE report the number of edges, of clusters and of taxa
  before and after agglomeration. Default to FALSE.

## Value

A
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object whose taxa are the network sequence clusters; or a data.frame
with columns `taxa` and `cluster` when `return_map = TRUE`; or an
`igraph` object when `return_graph = TRUE`.

## Details

vsearch is an external program and is not installed by phylopq. Use
[`MiscMetabar::is_vsearch_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_vsearch_installed.html)
to check for it and
[`MiscMetabar::install_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/reference/install_vsearch.html)
to install it.

`--allpairs_global` aligns every pair of sequences, so the run is
quadratic in the number of taxa. It is comfortable for a few thousand
taxa and slow well beyond that; subset or pre-cluster first in that
case.

Taxa that vsearch reports no hit for form their own single-taxon cluster
rather than being dropped, so no abundance is lost by the agglomeration.

## References

Forster D., Filker S., Kochems R., Breiner H.-W., Cordier T., Pawlowski
J., Stoeck T. (2019) A comparison of different ciliate metabarcode genes
as bioindicators for environmental impact assessments of salmon
aquaculture. *Environmental Microbiology* 21(11):4429-4443.
[doi:10.1111/1462-2920.14764](https://doi.org/10.1111/1462-2920.14764)

## See also

[`ssn_glom_scan_pq()`](https://adrientaudiere.github.io/phylopq/reference/ssn_glom_scan_pq.md),
[`phylo_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_pq.md),
[`MiscMetabar::postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.html),
[`MiscMetabar::merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.html)

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
phyloseq::ntaxa(df)
#> [1] 45

pq_nsc <- ssn_glom_pq(df, id = 0.9, verbose = TRUE)
#> ℹ 36 edges at id = 0.9 link 45 taxa into 28 clusters.
#> ✔ 45 taxa agglomerated into 28 network sequence clusters at id = 0.9.
phyloseq::ntaxa(pq_nsc)
#> [1] 28

# Inspect the taxa-to-cluster mapping without merging
map <- ssn_glom_pq(df, id = 0.9, return_map = TRUE)
head(map)
#>    taxa cluster
#> 1  ASV7       1
#> 2  ASV8       1
#> 3 ASV12       2
#> 4 ASV18       1
#> 5 ASV25       3
#> 6 ASV26       1

if (FALSE) { # \dontrun{
# The network itself, to plot or to analyse further
net <- ssn_glom_pq(df, id = 0.9, return_graph = TRUE)
igraph::V(net)$nsc
plot(net, vertex.label = NA, vertex.size = 5)

# Keep the taxonomy strictly conservative across merged taxa
pq_nsc <- ssn_glom_pq(df, id = 0.9, tax_adjust = 2)

# An explicit path to the executable
pq_nsc <- ssn_glom_pq(df, id = 0.9, vsearchpath = "/usr/local/bin/vsearch")
} # }
```
