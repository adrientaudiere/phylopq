# Assemble a reference package for phylogenetic placement

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Build the pair
[`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md)
needs — a reference **alignment** and a reference **tree** sharing the
same labels — from reference sequences that normally come from outside
the phyloseq object.

A reference package is the scientific heart of a placement: queries can
only be resolved as finely as the reference around the branch they land
on, so the reference is where the curation effort belongs. See the
*Choosing a reference package* section below.

## Usage

``` r
ref_package_pq(
  x,
  physeq = NULL,
  ref_taxa = NULL,
  tree = NULL,
  tree_method = c("nj", "upgma", "ml"),
  model = "GTR+G",
  align_method = c("decipher", "mafft"),
  mafft_exec = NULL,
  force_align = FALSE,
  drop_tips = TRUE,
  drop_seqs = FALSE,
  verbose = FALSE
)
```

## Arguments

- x:

  (required) The reference sequences: a path to a FASTA file, a
  `DNAStringSet`, or a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object whose `refseq` slot is used. Already-aligned sequences are kept
  as they are unless `force_align = TRUE`.

- physeq:

  Optional
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object to draw **additional** reference sequences from, on top of `x`.
  Combined with `ref_taxa`, this covers the case where a few
  well-identified taxa of your own dataset strengthen an external
  reference before the remaining taxa are placed on it.

- ref_taxa:

  Character vector of taxa names of `physeq` to add to the reference.
  Default to NULL, i.e. every taxon of `physeq`. Ignored when `physeq`
  is NULL.

- tree:

  Optional reference tree, a `phylo` object or a path to a Newick or
  Nexus file. **Give it whenever you have one**: a tree inferred from
  several loci, or constrained by morphology and ecology, carries far
  more evidence than one built here from the barcode alone. When NULL
  (default) a tree is inferred from the alignment with `tree_method`.

- tree_method:

  Method used when `tree` is NULL: `"nj"` (default, neighbour-joining on
  a maximum-likelihood distance), `"upgma"`, or `"ml"` (a
  maximum-likelihood tree optimised from the NJ starting tree, much
  slower). All three go through the `phangorn` package.

- model:

  The evolutionary model recorded in the package and handed to EPA-ng by
  [`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md),
  either a model string such as `"GTR+G"` (default) or the path to a
  RAxML-ng `.bestModel` file.

- align_method:

  Aligner used to build the alignment, `"decipher"` (default, pure R) or
  `"mafft"`. See
  [`MiscMetabar::align_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/align_pq.html).

- mafft_exec:

  Path to the MAFFT executable. Only used when `align_method = "mafft"`.

- force_align:

  Logical, if TRUE sequences that already share a width are realigned
  anyway. Default to FALSE.

- drop_tips:

  Logical, if TRUE (default) tips of `tree` with no sequence in the
  alignment are dropped. If FALSE their presence is an error.

- drop_seqs:

  Logical, if TRUE sequences with no tip in `tree` are dropped from the
  alignment. Default to FALSE, in which case they are an error — a
  sequence missing from the tree is usually a labelling mistake, not
  something to silently discard.

- verbose:

  Logical, if TRUE report the size of the package and every
  reconciliation step. Default to FALSE.

## Value

An object of class `phylopq_refpkg`, a list with the elements
`alignment` (a `DNAStringSet`), `tree` (a `phylo`), `model` and `n_ref`.
Pass it to
[`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md)
as `refpkg`.

## Choosing a reference package

EPA-ng places a query by asking where on the reference tree its sequence
fits best, so the alignment and the tree play two different roles and
answer to two different constraints.

- **The reference alignment must be the same marker as your queries** —
  the very locus your primers amplify, over a comparable region. Placing
  ITS2 reads on an ITS1 reference, or a short amplicon on a reference
  trimmed elsewhere, produces confident-looking nonsense. This is the
  constraint that cannot be relaxed.

- **The reference tree does not have to come from that marker.** It is
  taken as fixed and never re-estimated, so it is the right place to
  inject evidence the barcode does not carry: a multi-locus phylogeny, a
  topology constrained by morphology, ecology or mating tests, the
  published tree of an integrative-taxonomy study. A single short
  barcode rarely resolves deep relationships; a tree built from several
  loci and checked against morphological characters does, and placement
  inherits that resolution.

The usual arrangement is therefore an **external** reference package:
tips named for well-identified vouchers, the alignment restricted to
your marker, the topology taken from the best available phylogeny of the
group. Pass that tree through `tree` and only the alignment is built
here.

Tip labels and sequence names must match; use `drop_tips` and
`drop_seqs` to say what to do with the ones that do not.

## Three workflows

1.  **Fully external reference.** Sequences and tree both come from
    outside the phyloseq object; its taxa are the queries. This is the
    common case.

2.  **External sequences, tree inferred here.** No published tree is
    available for the group, so one is built from the reference
    sequences with `tree_method`. Weaker, but honest as long as the
    reference sequences are well identified.

3.  **Mixed reference.** A few securely identified taxa of your own
    dataset are added to an external reference with `physeq` and
    `ref_taxa`, to thicken the part of the tree your queries fall into;
    the remaining taxa are then placed with
    `place_pq(exclude_taxa = ref_taxa)`.

## See also

[`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md),
[`assign_placement_pq()`](https://adrientaudiere.github.io/phylopq/reference/assign_placement_pq.md),
[`MiscMetabar::align_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/align_pq.html),
[`add_tree_pq()`](https://adrientaudiere.github.io/phylopq/reference/add_tree_pq.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
library(MiscMetabar)
data(data_fungi_mini)
df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 9000)
#> Cleaning suppress 0 taxa (  ) and 6 sample(s) ( AD26-005-H_S10_MERGED.fastq.gz / CB8-019-H_S70_MERGED.fastq.gz / DY5-004-H_S97_MERGED.fastq.gz / N23-002-B_S130_MERGED.fastq.gz / NVABM0244-M_S137_MERGED.fastq.gz / T28-ABM602-B_S162_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 45
#> Number of filtered-out ASV 23
#> Number of kept ASV 22
#> Number of kept samples 131

# Workflow 2: reference sequences without a published tree.
# Here they are taken from the same object for the sake of the example;
# in practice they come from a curated database.
ref_taxa <- phyloseq::taxa_names(df)[1:8]
refpkg <- ref_package_pq(
  phyloseq::prune_taxa(ref_taxa, df),
  verbose = TRUE
)
#> Determining distance matrix based on shared 8-mers:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Clustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Aligning Sequences:
#> ================================================================================
#> 
#> Time difference of 0.04 secs
#> 
#> Iteration 1 of 2:
#> 
#> Determining distance matrix based on alignment:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Reclustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Realigning Sequences:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Alignment converged - skipping remaining iteration.
#> 
#> ✔ 8 sequences aligned with decipher over 363 positions.
#> ℹ Inferring a nj tree from the reference alignment; give `tree` to use a better
#>   one.
#> ✔ Reference package of 8 sequences over 363 positions.
refpkg
#> <phylopq reference package>
#>   8 reference sequences aligned over 363 positions
#>   a tree of 8 tips, with branch lengths
#>   model: GTR+G
# }

if (FALSE) { # \dontrun{
# Workflow 1: a fully external reference package, the common case.
# The tree comes from a multi-locus integrative study, the alignment from
# the same marker as your primers.
refpkg <- ref_package_pq(
  "unite_refs_ITS2.fasta",
  tree = "integrative_multilocus.newick",
  model = "reference.bestModel"
)
placements <- place_pq(data_fungi_mini, refpkg = refpkg)

# Workflow 3: strengthen an external reference with a few taxa of your own,
# then place the rest on it.
sure_taxa <- c("ASV1", "ASV12")
refpkg <- ref_package_pq(
  "unite_refs_ITS2.fasta",
  physeq = data_fungi_mini,
  ref_taxa = sure_taxa,
  align_method = "mafft"
)
placements <- place_pq(
  data_fungi_mini,
  refpkg = refpkg,
  exclude_taxa = sure_taxa
)

# A slower but better tree when no published one exists
refpkg <- ref_package_pq("refs.fasta", tree_method = "ml")
} # }
```
