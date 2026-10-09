# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository — the **phylopq** sub-package of the
pqverse.

## Package Overview

**phylopq** is the phylogenetic analysis layer of the pqverse. It
operates on `phyloseq` objects and provides tree construction from
taxonomy tables, phylogenetic distance metrics, phylogeny-aware
ordinations, and integration with external phylogenetic placement tools.

**Scope guard** (from `ROADMAP.md`):

- ✅ In scope: phylogeny-aware analyses on a single phyloseq object
  (tree from taxonomy, weighted/unweighted UniFrac-style distances,
  phylogeny-aware ordinations, tree-aware visualisations, integration
  with placement tools like epa-ng / BoSSA / gappa).
- ❌ Out of scope: pure ggplot2 wrappers (→ `ggplotpq`), data-structure
  utilities (→ `tidypq`), multi-phyloseq comparators (→ `comparpq`),
  general ML / networks / DAGs (→ `netaipq`), bootstrapping (→
  `bootpq`), reference DB I/O (→ `dbpq`).

**Dependency rule.** New heavy phylogenetic dependencies (e.g. `ape`
extensions, placement tools) live here. New pure-ggplot2 deps belong in
`ggplotpq`; new general-analysis deps belong in `netaipq`.

## Common Commands

``` bash
# Run code with loaded package
Rscript -e "devtools::load_all(); code"

# Run all tests
Rscript -e "devtools::test()"

# Run tests for files starting with {name}
Rscript -e "devtools::test(filter = '^{name}')"

# Generate documentation
Rscript -e "devtools::document()"

# Full package check
Rscript -e "devtools::check()"
```

## Coding Conventions

- Use base pipe (`|>`) not magrittr (`%>%`)
- Use `function() {}` for anonymous functions (not `\\()` for
  multi-statement)
- Line length limit: 120 characters
- Tests for `R/{name}.R` go in `tests/testthat/test_{name}.R`
  (underscore)
- Every user-facing function must be exported with full roxygen2
  documentation (`@param`, `@return`, `@export`, `@examples`, `@author`)
- Wrap roxygen comments at 80 characters
- CRAN example constraints: primary example in `\\donttest{}`, variants
  in `\\dontrun{}`; cap per-sample work at 5 samples via
  `prune_samples(sample_names(data_fungi_mini)[1:5], data_fungi_mini)`
- Guard every Suggests-package call with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) +
  [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
- Air format the package: `air format .` (then scope the diff — revert
  incidental reformats to unrelated files)

## Shipped so far

- [`taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.md)
  — build a `phylo` tree from a `tax_table()` (0.1.0, relocated from
  `comparpq`).

- [`add_tree_pq()`](https://adrientaudiere.github.io/phylopq/reference/add_tree_pq.md)
  — attach a tree (`phylo`, Newick/Nexus file, or
  [`taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.md)
  output) to a phyloseq object, reconciling tip labels and taxa names.

- [`phylo_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_pq.md)
  /
  [`phylo_glom_scan_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_scan_pq.md)
  — agglomerate taxa below a cophenetic distance threshold.

- [`delim_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_pq.md)
  /
  [`is_delim_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_delim_installed.md)
  — ABGD / ASAP species delimitation through `delimtools`.

- [`install_asap()`](https://adrientaudiere.github.io/phylopq/reference/install_asap.md)
  /
  [`install_abgd()`](https://adrientaudiere.github.io/phylopq/reference/install_abgd.md)
  — build the external ABGD / ASAP executables into
  `tools::R_user_dir("phylopq", "data")/bin`, where
  [`is_delim_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_delim_installed.md)
  now looks (option → user data dir → `PATH`). Upstream
  `bioinfo.mnhn.fr` is down, so sources come from the Internet Archive
  (ABGD, pristine) and the iTaxoTools mirrors (ASAP, whose `wrapio.h`
  Python shim is neutralised in a throwaway copy before `make`); `src`
  takes a URL, archive or directory for a manual install.

- [`align_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/align_pq.html)
  /
  [`is_mafft_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_mafft_installed.html)
  — **moved to MiscMetabar** and re-exported here (`R/align_pq.R` is now
  a re-export shim). MiscMetabar cannot depend on phylopq (golden rule)
  but everything depends on MiscMetabar, so this is what makes the MAFFT
  backend reachable from
  [`MiscMetabar::build_phytree_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/build_phytree_pq.html)
  and
  [`taxinfo::intra_taxnames_dist()`](https://adrientaudiere.github.io/taxinfo/reference/intra_taxnames_dist.html),
  which both gained `align_method`. The option is now
  `MiscMetabar.mafftpath`, not `phylopq.mafftpath`. Both call sites pass
  `force = TRUE` so an equal-width `refseq` is still aligned, as it was
  when they called DECIPHER directly.

- [`delim_multi_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_multi_pq.md)
  — scan ABGD across `slopes`, add one ASAP run, join the partitions and
  draw them along a tree with
  [`delimtools::delim_autoplot()`](https://legallab.github.io/delimtools/reference/delim_autoplot.html).
  Two upstream quirks are worked around: `delim_join()` strips every
  digit from delimitation names (so runs are submitted under digit-free
  aliases and renamed afterwards), and `delim_autoplot()` reads
  `posterior` / `support` node columns that a plain `phylo` lacks (so
  the tree is converted to `treedata` with both columns filled).

- [`ssn_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/ssn_glom_pq.md)
  /
  [`ssn_glom_scan_pq()`](https://adrientaudiere.github.io/phylopq/reference/ssn_glom_scan_pq.md)
  — network sequence clusters (NSC) from a `vsearch --allpairs_global`
  sequence-similarity network (Forster et al. 2019). Connected
  components are *single-linkage*, so a cluster can be far wider than
  `id`; the scan runs vsearch once at the lowest threshold and filters
  the identities at each of the others. `--strand` is **not** a valid
  `--allpairs_global` option; `--notrunclabels` is passed so taxa names
  survive the FASTA header.

- [`ref_package_pq()`](https://adrientaudiere.github.io/phylopq/reference/ref_package_pq.md)
  — assembles the reference alignment + tree that
  [`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md)
  needs. **The scientific point** (from the developer): the reference is
  normally *external* to the phyloseq object, and the two halves answer
  to different constraints — the alignment must be the same marker the
  primers amplify, whereas the tree is fixed and never re-estimated, so
  it should carry evidence the barcode cannot (multi- locus,
  morphology/ecology, integrative-taxonomy topologies). Three workflows
  are supported: fully external; external sequences with the tree
  inferred here (`tree_method`, weakest); mixed, where a few sure taxa
  of the dataset join the reference via `physeq`/`ref_taxa` and are then
  kept out of the queries with `place_pq(exclude_taxa = )`.
  Reconciliation is explicit: `drop_tips` yes by default, `drop_seqs`
  **no** by default (a sequence missing from the tree is a labelling
  bug, not something to discard), and an unrelated label set is named.

- [`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md)
  /
  [`assign_placement_pq()`](https://adrientaudiere.github.io/phylopq/reference/assign_placement_pq.md)
  — EPA-ng placement on a fixed reference tree, then
  `gappa examine assign` for the taxonomy. Takes either `refpkg` or
  `ref_alignment` + `ref_tree`, and `query_taxa` / `exclude_taxa` choose
  what gets placed. Queries are aligned into the reference alignment
  with `mafft --add --keeplength` (invoked directly:
  [`ips::mafft()`](https://rdrr.io/pkg/ips/man/mafft.html) offers no
  `--add`), and the `.jplace` is read with
  [`BoSSA::read_jplace()`](https://rdrr.io/pkg/BoSSA/man/read_jplace.html),
  which renumbers branches — hence both `edge` and `jplace_edge` in the
  result. `cmd_is_run = FALSE` builds the commands without the programs
  installed and without requiring the reference files to exist. Two
  gappa facts verified against v0.9.0: `per_query.tsv` is written
  **only** with `--per-query-results` (always passed; the default
  `profile.tsv` aggregates over the sample and has no query names), and
  it holds **one row per query per taxonomic depth** with `LWR`/`aLWR`
  columns — so the assignment is the deepest row clearing `min_alwr`
  (default 0.5), not the single row a naive reader would expect.
  Verified end-to-end against epa-ng v0.3.8 + gappa v0.9.0 with a
  leave-out reference package built from `data_fungi_mini` (see
  `test_place_pq.R`): 3 of 4 left-out queries were recovered exactly to
  genus, the fourth correctly to order.

- [`install_epang()`](https://adrientaudiere.github.io/phylopq/reference/install_epang.md)
  /
  [`install_gappa()`](https://adrientaudiere.github.io/phylopq/reference/install_gappa.md)
  — recursive `git clone` + cmake, because both carry their libraries as
  submodules that a release tarball lacks. Unlike ABGD/ASAP these cannot
  come from an archive.
  **[`install_epang()`](https://adrientaudiere.github.io/phylopq/reference/install_epang.md)
  additionally needs `bison` and `flex`**: epa-ng bundles `pll-modules`
  → `libpll`, whose CMake calls `BISON_TARGET()`, and without them the
  build dies at *configure* time with
  `Unknown CMake command "BISON_TARGET"` — several minutes into the
  clone, hence the up-front check in `check_place_build_tools()`. gappa
  builds on `genesis` alone and needs neither.

Demos: `arround_MiscMetabar/phylopq_demo.qmd` (0.2.0 features) and
`arround_MiscMetabar/phylopq_placement_demo.qmd` (NSC + placement).

## Remaining ROADMAP items

All three items that stood in the phylopq section of `ROADMAP.md` have
shipped. The section was refilled with three new ones, none of them
critical:

0.  **Vignette on the three
    [`ref_package_pq()`](https://adrientaudiere.github.io/phylopq/reference/ref_package_pq.md)
    workflows** — \[High/moderate\], and the one with a written plan
    already: `roadmap/place_pq_vignette.md` at the workspace root.
    Blocked on sourcing a real reference package (clade, published
    multi-locus tree, licensing); candidate clades and the sourcing gate
    are in that file. Every current example uses a synthetic or
    leave-out reference, which is exactly what the vignette must not do.
1.  Phylogenetic diversity metrics (Faith’s PD, MPD, MNTD, UniFrac) —
    \[Medium/moderate\]. `DESCRIPTION` advertises them and nothing in
    the pqverse provides them yet.
2.  Phylogeny-aware visualisation, a `plot_tree_pq()` on `ggtree` —
    \[Medium/moderate\].
3.  Placement follow-ups: `plot_placement_pq()` over BoSSA’s plotting,
    and a [`BoSSA::refpkg()`](https://rdrr.io/pkg/BoSSA/man/refpkg.html)
    reader so one `.refpkg` argument replaces `ref_alignment` +
    `ref_tree` + `model` — \[Low/moderate\].

## Cross-references

- Workspace CLAUDE.md: `pqverse/CLAUDE.md` (overall context)
- ROADMAP section:
  <https://github.com/adrientaudiere/pqverse/ROADMAP.md#phylopq--phylogenetic-analysis-for-phyloseq>
- Sister packages: `pqverse_pkg/MiscMetabar/`, `pqverse_pkg/bootpq/`,
  `pqverse_pkg/ggplotpq/`, `pqverse_pkg/tidypq/`, `pqverse_pkg/dbpq/`,
  `pqverse_pkg/comparpq/`, `pqverse_pkg/netaipq/`

## Agent skills

### Issue tracker

Issues and PRDs are tracked as GitHub issues via the `gh` CLI; external
PRs are not a triage surface. See `docs/agents/issue-tracker.md`.

### Triage labels

Uses the five canonical triage labels (`needs-triage`, `needs-info`,
`ready-for-agent`, `ready-for-human`, `wontfix`). See
`docs/agents/triage-labels.md`.

### Domain docs

Single-context: one `CONTEXT.md` + `docs/adr/` at the repo root. See
`docs/agents/domain.md`.
