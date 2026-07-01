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

## First migrations (priority order, from ROADMAP.md)

The current ROADMAP phylopq section is thin (only 2 items, one already
flagged as a `tidypq` candidate). The section likely needs enrichment
**before** the first feature batch — candidates to consider moving from
`MiscMetabar` (line 121, epa-ng / BoSSA / gappa placement) and `tidypq`
(line 161, phylo from taxonomy) once the home is settled.

1.  `taxonomy_to_phylo()` (working name) — build a `phylo` tree from a
    `tax_table()` — \[Medium/moderate\]. *Currently flagged as tidypq
    candidate in ROADMAP.md line 297.*
2.  ABGD / ASAP post-clustering wrapper — \[High/moderate\]. Prototype
    in `arround_MiscMetabar/delimtools_trying.R` (not yet present on
    disk).

See the R Feature Batch skill (`/r-feature-batch`) for the per-feature
workflow.

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
