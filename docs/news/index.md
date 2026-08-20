# Changelog

## phylopq 0.2.0 (Development version)

- Fix missing `Remotes` field in `DESCRIPTION` so that
  [`pak::pkg_install()`](https://pak.r-lib.org/reference/pkg_install.html)
  can resolve the GitHub-only dependency `MiscMetabar` when installing
  phylopq standalone.
- [`add_tree_pq()`](https://adrientaudiere.github.io/phylopq/reference/add_tree_pq.md)
  attaches a phylogenetic tree to a phyloseq object from a `phylo`
  object, a Newick or Nexus file, or directly from the taxonomy with
  [`taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.md),
  reconciling tip labels and taxa names explicitly (drop tips, prune
  taxa) and optionally rooting, ladderizing and computing arbitrary
  branch lengths.
- [`delim_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_pq.md)
  delimits species from the `refseq` slot with ABGD or ASAP through the
  `delimtools` package, aligning the sequences when needed and merging
  taxa of the same partition;
  [`is_delim_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_delim_installed.md)
  reports whether the external executable is available.
- [`phylo_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_pq.md)
  agglomerates taxa closer than a cophenetic distance threshold,
  returning either the merged phyloseq object or the taxa-to-cluster
  mapping, and
  [`phylo_glom_scan_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_scan_pq.md)
  reports the number of resulting taxa across a range of thresholds.

## phylopq 0.1.0

- [`taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.md)
  builds a “phylogenetic” `phylo` tree from the taxonomic ranks of a
  phyloseq object, with options to keep taxa names as tips
  (`use_taxa_names`) or collapse identical taxonomy paths, and to keep
  or drop singleton internal nodes (`internal_node_singletons`).
  Relocated from `comparpq`, its natural home in the pqverse.

## phylopq 0.0.0

- Initial development version of the package, providing phylogenetic
  analysis helpers for ‘phyloseq’ objects. Part of the ‘pqverse’
  ecosystem and built on top of ‘phyloseq’ and ‘ape’.
