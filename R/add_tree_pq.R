#' Attach, repair and prune a phylogenetic tree onto a phyloseq object
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Safely add a phylogenetic tree to the `phy_tree` slot of a
#' \code{\link[phyloseq]{phyloseq-class}} object. By default the tree is
#' supplied by the user as a `phylo` object (`ape` package) or as a path to a
#' Newick or Nexus file. Set `use_taxo_to_build_tree = TRUE` to build a
#' taxonomy-based tree on the fly with [taxo2tree()] instead.
#'
#' Tip labels and taxa names rarely match exactly. `add_tree_pq()` reconciles
#' them explicitly: tips absent from the phyloseq object are dropped, taxa
#' absent from the tree are pruned, and the function reports what it removed
#' instead of failing silently or letting `phyloseq` intersect quietly.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#' @param tree (required unless `use_taxo_to_build_tree` is TRUE) The tree to
#'   attach. Either a `phylo` object (`ape` package) or a length-one character
#'   giving the path to a Newick (`.nwk`, `.tre`, `.newick`) or Nexus (`.nex`,
#'   `.nexus`) file.
#' @param use_taxo_to_build_tree Logical, if TRUE a taxonomy-based tree is
#'   built from the `tax_table` slot with [taxo2tree()] and `tree` must be left
#'   NULL. Default to FALSE.
#' @param drop_tips Logical, if TRUE (default) tips of `tree` that are not
#'   taxa of `physeq` are dropped with [ape::drop.tip()].
#' @param prune_taxa Logical, if TRUE taxa of `physeq` that are not tips of
#'   `tree` are pruned with [phyloseq::prune_taxa()]. Default to FALSE, in
#'   which case such taxa trigger an error: dropping data is opt-in, since a
#'   tree that does not cover every taxon usually means something went wrong
#'   upstream.
#' @param root One of `"none"` (default, keep the tree as it is),
#'   `"midpoint"` (midpoint rooting, requires the `phangorn` package) or
#'   `"outgroup"` (root on `outgroup` with [ape::root()]).
#' @param outgroup Character vector of tip labels used as outgroup. Only used
#'   when `root = "outgroup"`.
#' @param ladderize Logical, if TRUE (default) the tree is ladderized with
#'   [ape::ladderize()].
#' @param compute_brlen Logical, if TRUE and the tree has no branch lengths,
#'   arbitrary branch lengths are computed with [ape::compute.brlen()]. This is
#'   required by downstream distance-based tools (e.g. [phylo_glom_pq()] or
#'   `phyloseq::UniFrac()`) but the resulting lengths carry no evolutionary
#'   meaning. Default to FALSE.
#' @param force Logical, if TRUE an existing `phy_tree` slot is replaced.
#'   Default to FALSE, in which case a non-empty slot triggers an error.
#' @param verbose Logical, if TRUE report the number of dropped tips and pruned
#'   taxa. Default to FALSE.
#' @param ... Additional arguments passed on to [taxo2tree()] when
#'   `use_taxo_to_build_tree` is TRUE.
#'
#' @return A \code{\link[phyloseq]{phyloseq-class}} object with a `phy_tree`
#'   slot whose tip labels match its taxa names exactly.
#'
#' @author Adrien Taudière
#'
#' @seealso [taxo2tree()], [phylo_glom_pq()],
#'   [MiscMetabar::build_phytree_pq()], [phyloseq::phy_tree()]
#'
#' @examples
#' \donttest{
#' library(MiscMetabar)
#' data(data_fungi_mini)
#'
#' # Attach a phylo object
#' tree <- taxo2tree(data_fungi_mini)
#' pq_tree <- add_tree_pq(data_fungi_mini, tree, verbose = TRUE)
#' pq_tree
#'
#' # Build a taxonomic tree and attach it in one call, adding arbitrary
#' # branch lengths
#' pq_brlen <- add_tree_pq(
#'   data_fungi_mini,
#'   use_taxo_to_build_tree = TRUE,
#'   compute_brlen = TRUE
#' )
#' ape::is.ultrametric(phyloseq::phy_tree(pq_brlen))
#'
#' # A tree covering only part of the taxa errors unless pruning is opt-in
#' small_tree <- ape::drop.tip(tree, tree$tip.label[1:10])
#' try(add_tree_pq(data_fungi_mini, small_tree))
#' pq_small <- add_tree_pq(
#'   data_fungi_mini,
#'   small_tree,
#'   prune_taxa = TRUE,
#'   verbose = TRUE
#' )
#' phyloseq::ntaxa(pq_small)
#' }
#'
#' \dontrun{
#' # Read a tree from a Newick file written by an external program
#' pq <- add_tree_pq(data_fungi_mini, "my_tree.nwk", root = "midpoint")
#'
#' # Replace an existing tree
#' pq <- add_tree_pq(pq, tree, force = TRUE)
#'
#' # Attach a real phylogeny inferred from the refseq slot with
#' # MiscMetabar::build_phytree_pq(), which returns a list of trees. Such a
#' # tree already carries branch lengths, so `compute_brlen` is not needed and
#' # the result can feed distance-based tools directly.
#' set.seed(22)
#' df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 9000)
#' phytree <- MiscMetabar::build_phytree_pq(df)
#'
#' pq_ml <- add_tree_pq(df, phytree$ML$tree, root = "midpoint")
#' phylo_glom_pq(pq_ml, h = 0.05, verbose = TRUE)
#'
#' # The other trees of the list are attached the same way
#' pq_upgma <- add_tree_pq(df, phytree$UPGMA)
#' pq_nj <- add_tree_pq(df, phytree$NJ, root = "midpoint")
#' }
#' @importFrom MiscMetabar verify_pq
#' @export
add_tree_pq <- function(
  physeq,
  tree = NULL,
  use_taxo_to_build_tree = FALSE,
  drop_tips = TRUE,
  prune_taxa = FALSE,
  root = c("none", "midpoint", "outgroup"),
  outgroup = NULL,
  ladderize = TRUE,
  compute_brlen = FALSE,
  force = FALSE,
  verbose = FALSE,
  ...
) {
  verify_pq(physeq)
  root <- match.arg(root)

  if (!is.null(physeq@phy_tree) && !force) {
    cli::cli_abort(c(
      "The {.arg physeq} object already contains a {.field phy_tree} slot.",
      "i" = "Use {.code force = TRUE} to replace it."
    ))
  }

  tree <- resolve_tree(physeq, tree, use_taxo_to_build_tree, ...)

  taxa <- phyloseq::taxa_names(physeq)
  shared <- intersect(tree$tip.label, taxa)

  if (length(shared) == 0) {
    first_tips <- utils::head(tree$tip.label, 3)
    first_taxa <- utils::head(taxa, 3)
    cli::cli_abort(c(
      "No tip label of {.arg tree} matches a taxa name of {.arg physeq}.",
      "i" = "First tips: {.val {first_tips}}.",
      "i" = "First taxa: {.val {first_taxa}}."
    ))
  }

  tips_to_drop <- setdiff(tree$tip.label, taxa)
  taxa_to_prune <- setdiff(taxa, tree$tip.label)

  if (length(tips_to_drop) > 0) {
    if (!drop_tips) {
      cli::cli_abort(c(
        "{length(tips_to_drop)} tip{?s} of {.arg tree} {?is/are} absent from
         {.arg physeq}.",
        "i" = "Use {.code drop_tips = TRUE} to drop them."
      ))
    }
    tree <- ape::drop.tip(tree, tips_to_drop)
  }

  if (length(taxa_to_prune) > 0) {
    if (!prune_taxa) {
      cli::cli_abort(c(
        "{length(taxa_to_prune)} taxa of {.arg physeq} {?is/are} absent from
         {.arg tree}.",
        "i" = "Use {.code prune_taxa = TRUE} to prune them."
      ))
    }
    physeq <- phyloseq::prune_taxa(shared, physeq)
  }

  if (compute_brlen && !has_brlen(tree)) {
    tree <- ape::compute.brlen(tree)
  }

  tree <- root_tree(tree, root, outgroup)

  if (ladderize) {
    tree <- ape::ladderize(tree)
  }

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Attached a tree with {length(tree$tip.label)} tip{?s}.",
      "i" = "{length(tips_to_drop)} tip{?s} dropped from the tree.",
      "i" = "{length(taxa_to_prune)} taxa pruned from the phyloseq object."
    ))
  }

  # Assign through phyloseq so every slot is reindexed on the tip order
  phyloseq::phy_tree(physeq) <- tree

  return(physeq)
}

#' Does a tree carry usable branch lengths?
#'
#' @description `phyloseq` turns a NULL `edge.length` into `numeric(0)`, so a
#'   plain `is.null()` check is not enough. Reaching `ape::dist.nodes()`
#'   without branch lengths segfaults.
#'
#' @param tree A `phylo` object.
#' @return A logical of length one.
#' @noRd
#' @keywords internal
has_brlen <- function(tree) {
  length(tree$edge.length) > 0
}

#' Resolve the `tree` argument of [add_tree_pq()]
#'
#' @inheritParams add_tree_pq
#' @param ... Additional arguments passed on to [taxo2tree()].
#' @return A `phylo` object.
#' @noRd
#' @keywords internal
resolve_tree <- function(physeq, tree, use_taxo_to_build_tree, ...) {
  if (use_taxo_to_build_tree) {
    if (!is.null(tree)) {
      cli::cli_abort(c(
        "{.arg tree} must be NULL when {.code use_taxo_to_build_tree = TRUE}.",
        "i" = "The tree is built from the {.field tax_table} slot."
      ))
    }
    return(taxo2tree(physeq, ...))
  }

  if (is.null(tree)) {
    cli::cli_abort(c(
      "{.arg tree} is required.",
      "i" = "Supply a {.cls phylo} object or the path to a tree file.",
      "i" = "Use {.code use_taxo_to_build_tree = TRUE} to build one from the
             {.field tax_table} slot with {.fn taxo2tree}."
    ))
  }

  if (inherits(tree, "phylo")) {
    return(tree)
  }

  if (methods::is(tree, "character") && length(tree) == 1) {
    if (!file.exists(tree)) {
      cli::cli_abort("The file {.file {tree}} does not exist.")
    }
    is_nexus <- grepl("\\.(nex|nexus)$", tree, ignore.case = TRUE)
    tree <- if (is_nexus) {
      ape::read.nexus(tree)
    } else {
      ape::read.tree(tree)
    }
    if (inherits(tree, "multiPhylo")) {
      tree <- tree[[1]]
    }
    return(tree)
  }

  cli::cli_abort(
    "{.arg tree} must be a {.cls phylo} object or a file path, not
     {.cls {class(tree)}}."
  )
}

#' Root a tree for [add_tree_pq()]
#'
#' @inheritParams add_tree_pq
#' @param tree A `phylo` object.
#' @return A `phylo` object.
#' @noRd
#' @keywords internal
root_tree <- function(tree, root, outgroup) {
  if (root == "none") {
    return(tree)
  }

  if (root == "midpoint") {
    if (!requireNamespace("phangorn", quietly = TRUE)) {
      cli::cli_abort(
        "Package {.pkg phangorn} is required for {.code root = \"midpoint\"}."
      )
    }
    if (!has_brlen(tree)) {
      cli::cli_abort(c(
        "Midpoint rooting requires branch lengths.",
        "i" = "Use {.code compute_brlen = TRUE} or supply a tree with branch
               lengths."
      ))
    }
    return(phangorn::midpoint(tree, node.labels = "label"))
  }

  if (is.null(outgroup)) {
    cli::cli_abort(
      "{.arg outgroup} must be given when {.code root = \"outgroup\"}."
    )
  }
  missing_tips <- setdiff(outgroup, tree$tip.label)
  if (length(missing_tips) > 0) {
    cli::cli_abort(
      "{length(missing_tips)} {.arg outgroup} tip{?s} {?is/are} absent from the
       tree: {.val {missing_tips}}."
    )
  }
  ape::root(tree, outgroup = outgroup, resolve.root = TRUE)
}
