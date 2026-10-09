# Assign a taxonomy from a phylogenetic placement

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Turn the `.jplace` file produced by
[`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md)
into a taxonomy, with `gappa examine assign` (Czech et al. 2020). Each
query is assigned the consensus taxonomic path of the reference tips
around the branch it was placed on, so a query that falls deep inside a
well-sampled clade is resolved to a fine rank while one placed near the
root is resolved only to a coarse one.

This is the complement of a similarity-based assignment: it answers
"where does this sequence sit in the reference phylogeny", not "which
reference sequence does it look most like".

## Usage

``` r
assign_placement_pq(
  jplace,
  taxon_file,
  physeq = NULL,
  ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  min_alwr = 0.5,
  exec = NULL,
  outdir = NULL,
  gappa_args = "",
  cmd_is_run = TRUE,
  verbose = FALSE
)
```

## Arguments

- jplace:

  (required) Path to a `.jplace` file, or a data.frame returned by
  [`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md)
  (its `jplace_file` attribute is used).

- taxon_file:

  (required) Path to the tab-separated file mapping each reference tip
  label to its taxonomic path, in the format `gappa examine assign`
  expects: one line per tip, the label, a tab, then the path with its
  ranks separated by `;`.

- physeq:

  Optional
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object. When given, the assigned taxonomy is written into its
  `tax_table` and the object is returned instead of the table.

- ranks:

  Character vector naming the ranks of `taxon_file`, used as the columns
  of the returned table. Default to the seven usual ranks.

- min_alwr:

  Numeric in \[0, 1\], the accumulated likelihood weight a clade must
  reach for a query to be assigned that deep. Default to 0.5. gappa
  reports one row per query **and per taxonomic depth**, each with the
  weight accumulated over that clade; the assignment kept here is the
  deepest row still clearing this floor, so raising `min_alwr` yields
  coarser but safer assignments and lowering it yields finer but shakier
  ones.

- exec:

  Path to the gappa executable. Default to NULL, in which case it is
  looked up with
  [`is_gappa_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_gappa_installed.md):
  first the `phylopq.gappapath` option, then a copy installed by
  [`install_gappa()`](https://adrientaudiere.github.io/phylopq/reference/install_gappa.md),
  then the system `PATH`.

- outdir:

  Directory gappa writes into. Default to NULL, i.e. a temporary
  location.

- gappa_args:

  A length-one character vector of further arguments passed on to gappa,
  e.g. `"--consensus-thresh 0.9"`. Default to `""`.

- cmd_is_run:

  Logical, if FALSE the command is returned instead of being run.
  Default to TRUE.

- verbose:

  Logical, if TRUE report each step and let gappa write to the console.
  Default to FALSE.

## Value

A data.frame with one row per query and the columns `query`, `alwr` (the
accumulated likelihood weight of the rank it was assigned to) and one
per rank of `ranks`; or the `physeq` object with that taxonomy in its
`tax_table` when `physeq` is given; or the command as a length-one
character vector when `cmd_is_run = FALSE`.

## Details

gappa is an external program and is not installed with phylopq. The
quickest route is conda (`conda install -c bioconda gappa`);
[`install_gappa()`](https://adrientaudiere.github.io/phylopq/reference/install_gappa.md)
builds it from source instead.

Queries that gappa returns no assignment for keep `NA` at every rank
rather than being dropped, so the returned table always covers every
query of the `.jplace` file.

The `--per-query-results` flag is always passed, because gappa writes
`per_query.tsv` only when asked to and this function reads nothing else;
`profile.tsv`, which it writes by default, aggregates over the whole
sample and carries no query names.

## References

Czech L., Barbera P., Stamatakis A. (2020) Genesis and Gappa:
processing, analyzing and visualizing phylogenetic (placement) data.
*Bioinformatics* 36(10):3263-3265.
[doi:10.1093/bioinformatics/btaa070](https://doi.org/10.1093/bioinformatics/btaa070)

## See also

[`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md),
[`is_gappa_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_gappa_installed.md),
[`install_gappa()`](https://adrientaudiere.github.io/phylopq/reference/install_gappa.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
# The command is built without gappa being installed
assign_placement_pq(
  jplace = "epa_result.jplace",
  taxon_file = "reference_taxonomy.tsv",
  cmd_is_run = FALSE
)
#> [1] "'gappa' examine assign --jplace-path 'epa_result.jplace' --taxon-file 'reference_taxonomy.tsv' --out-dir '/tmp/RtmpSIJqWk/phylopq_assign_c350d7ed18da2' --per-query-results --allow-file-overwriting "
# }

if (FALSE) { # \dontrun{
taxo <- assign_placement_pq(
  jplace = "epa_result.jplace",
  taxon_file = "reference_taxonomy.tsv",
  verbose = TRUE
)
head(taxo)

# Only keep assignments the placement is confident about
taxo[!is.na(taxo$alwr) & taxo$alwr >= 0.9, ]

# Coarser but safer: a clade must carry 90% of the weight to be used
taxo_safe <- assign_placement_pq(
  jplace = "epa_result.jplace",
  taxon_file = "reference_taxonomy.tsv",
  min_alwr = 0.9
)

# Write the assignment straight into a phyloseq object
pq <- assign_placement_pq(
  jplace = placements,
  taxon_file = "reference_taxonomy.tsv",
  physeq = df
)
phyloseq::tax_table(pq)[1:5, ]

# A stricter consensus
taxo <- assign_placement_pq(
  jplace = "epa_result.jplace",
  taxon_file = "reference_taxonomy.tsv",
  gappa_args = "--consensus-thresh 0.9"
)
} # }
```
