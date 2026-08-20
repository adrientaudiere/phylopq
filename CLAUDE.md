# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository — the **phylopq** sub-package of the
pqverse.

## Package Overview

**phylopq** is the phylogenetic analysis layer of the pqverse. It
operates on `phyloseq` objects and provides tree construction from
taxonomy tables, phylogenetic distance metrics, phylogeny-aware
ordinations, and integration with external phylogenetic placement
tools.

**Scope guard** (from `ROADMAP.md`):

- ✅ In scope: phylogeny-aware analyses on a single phyloseq object
  (tree from taxonomy, weighted/unweighted UniFrac-style distances,
  phylogeny-aware ordinations, tree-aware visualisations, integration
  with placement tools like epa-ng / BoSSA / gappa).
- ❌ Out of scope: pure ggplot2 wrappers (→ `ggplotpq`), data-structure
  utilities (→ `tidypq`), multi-phyloseq comparators (→ `comparpq`),
  general ML / networks / DAGs (→ `netaipq`), bootstrapping
  (→ `bootpq`), reference DB I/O (→ `dbpq`).

**Dependency rule.** New heavy phylogenetic dependencies (e.g. `ape`
extensions, placement tools) live here. New pure-ggplot2 deps belong in
`ggplotpq`; new general-analysis deps belong in `netaipq`.

## Common Commands

```bash
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
- Use `function() {}` for anonymous functions (not `\\()` for multi-statement)
- Line length limit: 120 characters
- Tests for `R/{name}.R` go in `tests/testthat/test_{name}.R` (underscore)
- Every user-facing function must be exported with full roxygen2
  documentation (`@param`, `@return`, `@export`, `@examples`,
  `@author`)
- Wrap roxygen comments at 80 characters
- CRAN example constraints: primary example in `\\donttest{}`, variants in
  `\\dontrun{}`; cap per-sample work at 5 samples via
  `prune_samples(sample_names(data_fungi_mini)[1:5], data_fungi_mini)`
- Guard every Suggests-package call with `requireNamespace()` +
  `cli::cli_abort()`
- Air format the package: `air format .` (then scope the diff — revert
  incidental reformats to unrelated files)

## Shipped so far

- `taxo2tree()` — build a `phylo` tree from a `tax_table()` (0.1.0,
  relocated from `comparpq`).
- `add_tree_pq()` — attach a tree (`phylo`, Newick/Nexus file, or
  `taxo2tree()` output) to a phyloseq object, reconciling tip labels
  and taxa names.
- `phylo_glom_pq()` / `phylo_glom_scan_pq()` — agglomerate taxa below a
  cophenetic distance threshold.
- `delim_pq()` / `is_delim_installed()` — ABGD / ASAP species
  delimitation through `delimtools`.

Demo: `arround_MiscMetabar/phylopq_demo.qmd`.

## Remaining ROADMAP items

The phylopq section of `ROADMAP.md` is thin and holds no `easy` item.
Enrich it before the next feature batch. What is left:

1. Phylogenetic placement via epa-ng / BoSSA / gappa — [High/hard].
2. Sequence similarity networks / NSC reclustering — [Low/hard].

Note that the `(source:)` pointers of both are stale.

See the `/pqverse-add-features` skill for the per-feature workflow.

## Cross-references

- Workspace CLAUDE.md: `pqverse/CLAUDE.md` (overall context)
- ROADMAP section: <https://github.com/adrientaudiere/pqverse/ROADMAP.md#phylopq--phylogenetic-analysis-for-phyloseq>
- Sister packages: `pqverse_pkg/MiscMetabar/`, `pqverse_pkg/bootpq/`,
  `pqverse_pkg/ggplotpq/`, `pqverse_pkg/tidypq/`, `pqverse_pkg/dbpq/`,
  `pqverse_pkg/comparpq/`, `pqverse_pkg/netaipq/`

## Agent skills

### Issue tracker

Issues and PRDs are tracked as GitHub issues via the `gh` CLI; external PRs are not a triage surface. See `docs/agents/issue-tracker.md`.

### Triage labels

Uses the five canonical triage labels (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`). See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: one `CONTEXT.md` + `docs/adr/` at the repo root. See `docs/agents/domain.md`.
