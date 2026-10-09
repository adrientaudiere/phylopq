# Package index

## Trees

Build a tree from the taxonomy, attach it to a ‘phyloseq’ object and
agglomerate taxa along it.

- [`add_tree_pq()`](https://adrientaudiere.github.io/phylopq/reference/add_tree_pq.md)
  : Attach, repair and prune a phylogenetic tree onto a phyloseq object
- [`phylo_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_pq.md)
  : Agglomerate taxa closer than a cophenetic distance threshold
- [`phylo_glom_scan_pq()`](https://adrientaudiere.github.io/phylopq/reference/phylo_glom_scan_pq.md)
  : Scan several cophenetic thresholds before agglomerating
- [`taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.md)
  : Convert taxonomy dataframe to phylogenetic tree

## Species delimitation and clustering

Regroup taxa by a barcode gap (ABGD, ASAP) or by a sequence-similarity
network (NSC).

- [`delim_multi_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_multi_pq.md)
  : Compare several species delimitations of a phyloseq object
- [`delim_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_pq.md)
  : Species delimitation of a phyloseq object with ABGD or ASAP
- [`ssn_glom_pq()`](https://adrientaudiere.github.io/phylopq/reference/ssn_glom_pq.md)
  : Agglomerate taxa into network sequence clusters (NSC)
- [`ssn_glom_scan_pq()`](https://adrientaudiere.github.io/phylopq/reference/ssn_glom_scan_pq.md)
  : Scan several similarity thresholds before agglomerating

## Phylogenetic placement

Place reference sequences on a fixed reference phylogeny with EPA-ng and
turn the placement into a taxonomy with gappa.

- [`assign_placement_pq()`](https://adrientaudiere.github.io/phylopq/reference/assign_placement_pq.md)
  : Assign a taxonomy from a phylogenetic placement
- [`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md)
  : Place the reference sequences of a phyloseq object on a reference
  tree
- [`print(`*`<phylopq_refpkg>`*`)`](https://adrientaudiere.github.io/phylopq/reference/print.phylopq_refpkg.md)
  : Print a reference package
- [`ref_package_pq()`](https://adrientaudiere.github.io/phylopq/reference/ref_package_pq.md)
  : Assemble a reference package for phylogenetic placement

## External programs

Check for and install the command-line programs the functions above
drive.

- [`install_abgd()`](https://adrientaudiere.github.io/phylopq/reference/install_abgd.md)
  : Install the ABGD species-delimitation program
- [`install_asap()`](https://adrientaudiere.github.io/phylopq/reference/install_asap.md)
  : Install the ASAP species-delimitation program
- [`install_epang()`](https://adrientaudiere.github.io/phylopq/reference/install_epang.md)
  : Install the EPA-ng phylogenetic placement program
- [`install_gappa()`](https://adrientaudiere.github.io/phylopq/reference/install_gappa.md)
  : Install the gappa placement-analysis program
- [`is_delim_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_delim_installed.md)
  : Is an ABGD or ASAP executable available?
- [`is_epang_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_epang_installed.md)
  : Is EPA-ng available?
- [`is_gappa_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_gappa_installed.md)
  : Is gappa available?

## Re-exported from MiscMetabar

Alignment helpers that live in MiscMetabar so that every pqverse package
can reach them, and are re-exported here.

- [`reexports`](https://adrientaudiere.github.io/phylopq/reference/reexports.md)
  [`align_pq`](https://adrientaudiere.github.io/phylopq/reference/reexports.md)
  [`is_mafft_installed`](https://adrientaudiere.github.io/phylopq/reference/reexports.md)
  : Objects re-exported from MiscMetabar
