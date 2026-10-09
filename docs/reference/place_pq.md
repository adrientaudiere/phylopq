# Place the reference sequences of a phyloseq object on a reference tree

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Run a *phylogenetic placement* of the `refseq` slot of a
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object onto a fixed reference phylogeny with EPA-ng (Barbera et al.
2019). Rather than building a tree from the metabarcoding reads
themselves — which short, fast-evolving markers rarely support — each
query sequence is grafted onto the branch of a curated reference tree it
fits best, with a likelihood weight attached.

The three steps are chained for you:

1.  the query sequences are aligned **into** the reference alignment
    with `mafft --add --keeplength`, so that they end up on exactly the
    reference's columns, which EPA-ng requires;

2.  EPA-ng places them on the reference tree and writes a `.jplace`
    file;

3.  the `.jplace` file is read back with
    [`BoSSA::read_jplace()`](https://rdrr.io/pkg/BoSSA/man/read_jplace.html)
    and flattened into a data.frame.

Use
[`ref_package_pq()`](https://adrientaudiere.github.io/phylopq/reference/ref_package_pq.md)
to assemble the reference and
[`assign_placement_pq()`](https://adrientaudiere.github.io/phylopq/reference/assign_placement_pq.md)
to turn the placements into a taxonomy.

## Usage

``` r
place_pq(
  physeq,
  ref_alignment = NULL,
  ref_tree = NULL,
  refpkg = NULL,
  model = "GTR+G",
  query_taxa = NULL,
  exclude_taxa = NULL,
  query = NULL,
  exec = NULL,
  mafft_exec = NULL,
  outdir = NULL,
  threads = 1,
  epang_args = "",
  best_only = TRUE,
  return_jplace = FALSE,
  keep_temporary_files = FALSE,
  cmd_is_run = TRUE,
  verbose = FALSE
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object with a non-empty `refseq` slot, or a `DNAStringSet` of query
  sequences.

- ref_alignment:

  The reference alignment, either a path to an aligned FASTA file or a
  `DNAStringSet` whose sequences all share a width. It must cover the
  **same marker as your queries**. Required unless `refpkg` is given.

- ref_tree:

  The reference tree, either a path to a Newick file or a `phylo`
  object. Its tip labels must match the names of `ref_alignment`. It is
  taken as fixed and never re-estimated, so it may — and ideally should
  — come from more evidence than the barcode alone. Required unless
  `refpkg` is given.

- refpkg:

  A `phylopq_refpkg` object from
  [`ref_package_pq()`](https://adrientaudiere.github.io/phylopq/reference/ref_package_pq.md),
  carrying the reference alignment, tree and model together. Give either
  this or `ref_alignment` + `ref_tree`.

- model:

  The evolutionary model handed to EPA-ng, either a model string such as
  `"GTR+G"` (default) or the path to a RAxML-ng `.bestModel` file. A
  model fitted on the reference tree is strongly preferable to the
  default.

- query_taxa:

  Character vector of taxa of `physeq` to place. Default to NULL, i.e.
  all of them.

- exclude_taxa:

  Character vector of taxa of `physeq` **not** to place. Use it for the
  taxa you put into the reference package itself, which would otherwise
  be placed on a branch they define.

- query:

  Optional path to a FASTA file of query sequences **already aligned**
  on the reference's columns. Default to NULL, i.e. the `refseq` slot of
  `physeq` is aligned for you with MAFFT.

- exec:

  Path to the EPA-ng executable. Default to NULL, in which case it is
  looked up with
  [`is_epang_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_epang_installed.md):
  first the `phylopq.epangpath` option, then a copy installed by
  [`install_epang()`](https://adrientaudiere.github.io/phylopq/reference/install_epang.md),
  then the system `PATH`.

- mafft_exec:

  Path to the MAFFT executable, used to align the queries. Default to
  NULL, i.e. the `MiscMetabar.mafftpath` option then the system `PATH`.
  Ignored when `query` is given.

- outdir:

  Directory EPA-ng writes into. Default to NULL, i.e. a temporary
  location. The `.jplace` file is named `epa_result.jplace` there.

- threads:

  Number of threads given to EPA-ng and to MAFFT. Default to 1.

- epang_args:

  A length-one character vector of further arguments passed on to
  EPA-ng, e.g. `"--filter-acc-lwr 0.99"`. Default to `""`.

- best_only:

  Logical, if TRUE (default) only the placement with the highest
  likelihood weight ratio is kept for each query. If FALSE every
  candidate branch is returned.

- return_jplace:

  Logical, if TRUE return the `jplace` object of
  [`BoSSA::read_jplace()`](https://rdrr.io/pkg/BoSSA/man/read_jplace.html)
  itself instead of the flattened data.frame. Default to FALSE.

- keep_temporary_files:

  Logical, if TRUE the FASTA files handed to MAFFT and EPA-ng are kept.
  Default to FALSE.

- cmd_is_run:

  Logical, if FALSE the commands are returned as a named character
  vector instead of being run. Default to TRUE. Useful to inspect, log
  or adapt the invocation, and to document the workflow without the
  external programs installed.

- verbose:

  Logical, if TRUE report each step and let EPA-ng and MAFFT write to
  the console. Default to FALSE.

## Value

A data.frame with one row per placement and the columns `query` (the
taxa name), `edge` (the branch index in the returned tree),
`jplace_edge` (the branch number as EPA-ng numbered it),
`like_weight_ratio`, `likelihood`, `distal_length` and `pendant_length`,
carrying the path to the `.jplace` file as its `jplace_file` attribute
and the reference tree as its `tree` attribute.

A `jplace` object when `return_jplace = TRUE`, or a named character
vector of commands when `cmd_is_run = FALSE`.

## Details

EPA-ng is an external program and is not installed with phylopq. The
quickest route is conda (`conda install -c bioconda epa-ng`);
[`install_epang()`](https://adrientaudiere.github.io/phylopq/reference/install_epang.md)
builds it from source instead.

The likelihood weight ratio (LWR) of a placement is the share of the
total likelihood that branch accounts for; it is the number to filter
on. A query spread thinly over many branches is placed with little
confidence, whatever its best branch is, which is why
`best_only = FALSE` is worth inspecting before trusting an assignment.

## The reference is the science

Placement never questions the reference tree — it only asks where a
query fits on it. Everything therefore rests on that reference, and the
two halves of it answer to different constraints:

- **The reference alignment must match your primers.** It has to be the
  same locus, over a comparable region, as the amplicon you sequenced.
  This is the constraint that cannot be relaxed: placing ITS2 reads on
  an ITS1 reference yields confident-looking nonsense.

- **The reference tree need not come from that locus at all.** Being
  fixed, it is the natural place to inject evidence a single barcode
  cannot carry — a multi-locus phylogeny, a topology constrained by
  morphology, ecology or mating data, the published tree of an
  integrative-taxonomy study. Short markers rarely resolve deep
  relationships; a tree built from several loci and morpho-ecological
  characters does, and the placement inherits that resolution.

In practice the reference package is almost always **external** to the
phyloseq object: curated vouchers, not your own ASVs. Three arrangements
are supported, all through
[`ref_package_pq()`](https://adrientaudiere.github.io/phylopq/reference/ref_package_pq.md):

1.  a fully external reference — sequences *and* tree from a database or
    a published study, every taxon of `physeq` a query;

2.  external reference sequences with the tree inferred from them, when
    no published phylogeny exists for the group;

3.  a mixed reference — a few securely identified taxa of your own
    dataset added to an external reference to thicken the clade your
    queries fall into, the rest placed with `exclude_taxa`.

## References

Barbera P., Kozlov A.M., Czech L., Morel B., Darriba D., Flouri T.,
Stamatakis A. (2019) EPA-ng: Massively Parallel Evolutionary Placement
of Genetic Sequences. *Systematic Biology* 68(2):365-369.
[doi:10.1093/sysbio/syy054](https://doi.org/10.1093/sysbio/syy054)

Matsen F.A., Kodner R.B., Armbrust E.V. (2010) pplacer: linear time
maximum-likelihood and Bayesian phylogenetic placement of sequences onto
a fixed reference tree. *BMC Bioinformatics* 11:538.
[doi:10.1186/1471-2105-11-538](https://doi.org/10.1186/1471-2105-11-538)

## See also

[`ref_package_pq()`](https://adrientaudiere.github.io/phylopq/reference/ref_package_pq.md),
[`assign_placement_pq()`](https://adrientaudiere.github.io/phylopq/reference/assign_placement_pq.md),
[`is_epang_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_epang_installed.md),
[`install_epang()`](https://adrientaudiere.github.io/phylopq/reference/install_epang.md),
[`MiscMetabar::align_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/align_pq.html),
[`BoSSA::read_jplace()`](https://rdrr.io/pkg/BoSSA/man/read_jplace.html)

## Author

Adrien Taudière

## Examples

``` r
# The commands are built without EPA-ng or MAFFT being installed
# \donttest{
library(MiscMetabar)
data(data_fungi_mini)
df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 12000)
#> Cleaning suppress 0 taxa (  ) and 7 sample(s) ( AD26-005-H_S10_MERGED.fastq.gz / CB8-019-H_S70_MERGED.fastq.gz / DY5-004-H_S97_MERGED.fastq.gz / F7-015-M_S106_MERGED.fastq.gz / N23-002-B_S130_MERGED.fastq.gz / NVABM0244-M_S137_MERGED.fastq.gz / T28-ABM602-B_S162_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 45
#> Number of filtered-out ASV 32
#> Number of kept ASV 13
#> Number of kept samples 130

place_pq(
  df,
  ref_alignment = "reference_alignment.fasta",
  ref_tree = "reference_tree.newick",
  model = "GTR+G",
  cmd_is_run = FALSE
)
#>                                                                                                                                                                                                                                      mafft 
#>                                              "'mafft' --add '/tmp/RtmpSIJqWk/phylopq_query_raw_c350d4f19af71.fasta' --keeplength --thread 1 --quiet 'reference_alignment.fasta' > '/tmp/RtmpSIJqWk/phylopq_query_add_c350d63f04bff.fasta'" 
#>                                                                                                                                                                                                                                     epa-ng 
#> "'epa-ng' --ref-msa 'reference_alignment.fasta' --tree 'reference_tree.newick' --query '/tmp/RtmpSIJqWk/phylopq_query_ali_c350d5cd1ebf6.fasta' --model 'GTR+G' --outdir '/tmp/RtmpSIJqWk/phylopq_place_c350d2c175fe0' --threads 1 --redo " 
# }

if (FALSE) { # \dontrun{
# Workflow 1 - a fully external reference package: curated sequences for
# your marker, and a tree from a multi-locus integrative study.
refpkg <- ref_package_pq(
  "unite_refs_ITS2.fasta",
  tree = "integrative_multilocus.newick",
  model = "reference.bestModel"
)
placements <- place_pq(df, refpkg = refpkg, threads = 4, verbose = TRUE)
head(placements)
attr(placements, "jplace_file")

# Workflow 2 - external sequences, no published tree for the group
refpkg <- ref_package_pq("refs.fasta", tree_method = "ml")
placements <- place_pq(df, refpkg = refpkg)

# Workflow 3 - strengthen the reference with a few taxa you are sure of,
# then place only the remaining ones
sure_taxa <- c("ASV1", "ASV12")
refpkg <- ref_package_pq(
  "refs.fasta",
  physeq = df,
  ref_taxa = sure_taxa
)
placements <- place_pq(df, refpkg = refpkg, exclude_taxa = sure_taxa)

# The alignment and the tree can still be given separately
placements <- place_pq(
  df,
  ref_alignment = "reference_alignment.fasta",
  ref_tree = "reference_tree.newick"
)

# Every candidate branch, to judge how confident each placement is
all_places <- place_pq(
  df,
  ref_alignment = "reference_alignment.fasta",
  ref_tree = "reference_tree.newick",
  best_only = FALSE
)

# Discard the queries EPA-ng is not confident about
placements[placements$like_weight_ratio >= 0.9, ]

# The jplace object itself, for BoSSA's plotting functions
jp <- place_pq(
  df,
  ref_alignment = "reference_alignment.fasta",
  ref_tree = "reference_tree.newick",
  return_jplace = TRUE
)
plot(jp, type = "precise")
} # }
```
