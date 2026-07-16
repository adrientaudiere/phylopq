# phylopq 0.2.0 (Development version)

* Fix missing `Remotes` field in `DESCRIPTION` so that `pak::pkg_install()` can resolve the GitHub-only dependency `MiscMetabar` when installing phylopq standalone.

# phylopq 0.1.0
* `taxo2tree()` builds a "phylogenetic" `phylo` tree from the taxonomic ranks of a phyloseq object, with options to keep taxa names as tips (`use_taxa_names`) or collapse identical taxonomy paths, and to keep or drop singleton internal nodes (`internal_node_singletons`). Relocated from `comparpq`, its natural home in the pqverse.

# phylopq 0.0.0
* Initial development version of the package, providing phylogenetic
  analysis helpers for 'phyloseq' objects. Part of the 'pqverse'
  ecosystem and built on top of 'phyloseq' and 'ape'.
