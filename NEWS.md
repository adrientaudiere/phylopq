# phylopq 0.2.0 (Development version)

* Fix missing `Remotes` field in `DESCRIPTION` so that `pak::pkg_install()` can resolve the GitHub-only dependency `MiscMetabar` when installing phylopq standalone.
* `add_tree_pq()` attaches a phylogenetic tree to a phyloseq object from a `phylo` object or a Newick or Nexus file, or from the taxonomy with `taxo2tree()` when `use_taxo_to_build_tree = TRUE`, reconciling tip labels and taxa names explicitly (drop tips, prune taxa) and optionally rooting, ladderizing and computing arbitrary branch lengths.
* `delim_pq()` delimits species from the `refseq` slot with ABGD or ASAP through the `delimtools` package, aligning the sequences when needed and merging taxa of the same partition; `is_delim_installed()` reports whether the external executable is available.
* `phylo_glom_pq()` agglomerates taxa closer than a cophenetic distance threshold, returning either the merged phyloseq object or the taxa-to-cluster mapping, and `phylo_glom_scan_pq()` reports the number of resulting taxa across a range of thresholds.

# phylopq 0.1.0
* `taxo2tree()` builds a "phylogenetic" `phylo` tree from the taxonomic ranks of a phyloseq object, with options to keep taxa names as tips (`use_taxa_names`) or collapse identical taxonomy paths, and to keep or drop singleton internal nodes (`internal_node_singletons`). Relocated from `comparpq`, its natural home in the pqverse.

# phylopq 0.0.0
* Initial development version of the package, providing phylogenetic
  analysis helpers for 'phyloseq' objects. Part of the 'pqverse'
  ecosystem and built on top of 'phyloseq' and 'ape'.
