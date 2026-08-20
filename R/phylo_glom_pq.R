#' Agglomerate taxa closer than a cophenetic distance threshold
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Merge taxa (ASV, OTU, ...) whose patristic (cophenetic)
#' distance on the `phy_tree` of a
#' \code{\link[phyloseq]{phyloseq-class}} object is below a threshold `h`.
#' Taxa are clustered with [stats::hclust()] on the cophenetic distance
#' matrix, cut at height `h` with [stats::cutree()], and merged with
#' [MiscMetabar::merge_taxa_vec()].
#'
#' This is a faster and more transparent alternative to
#' [phyloseq::tip_glom()]: the clustering method is explicit, the mapping
#' from taxa to clusters can be returned instead of the merged object, and
#' the number of merged taxa is reported.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   with a non-empty `phy_tree` slot carrying branch lengths.
#' @param h (required) Numeric, the cophenetic distance threshold at which the
#'   tree-based dendrogram is cut. Taxa within `h` of each other are merged.
#' @param hclust_method Agglomeration method passed on to [stats::hclust()].
#'   Default to `"average"` (UPGMA), as in [phyloseq::tip_glom()].
#' @param tax_adjust Handling of taxonomic disagreements within a cluster. See
#'   [MiscMetabar::merge_taxa_vec()]. Default to 1 (phyloseq-compatible).
#' @param rank_propagation Logical, default TRUE, whether bad ranks are
#'   propagated to lower ranks. See [MiscMetabar::merge_taxa_vec()].
#' @param return_map Logical, if TRUE return a data.frame mapping each taxon to
#'   its cluster instead of the agglomerated phyloseq object. Default to FALSE.
#' @param verbose Logical, if TRUE report the number of taxa before and after
#'   agglomeration. Default to FALSE.
#'
#' @return A \code{\link[phyloseq]{phyloseq-class}} object with agglomerated
#'   taxa, or a data.frame with columns `taxa` and `cluster` when
#'   `return_map = TRUE`.
#'
#' @author Adrien Taudière
#'
#' @seealso [phylo_glom_scan_pq()], [add_tree_pq()],
#'   [MiscMetabar::merge_taxa_vec()], [phyloseq::tip_glom()]
#'
#' @examples
#' \donttest{
#' library(MiscMetabar)
#' data(data_fungi_mini)
#'
#' pq <- add_tree_pq(
#'   data_fungi_mini,
#'   use_taxo_to_build_tree = TRUE,
#'   compute_brlen = TRUE
#' )
#' phyloseq::ntaxa(pq)
#'
#' pq_glom <- phylo_glom_pq(pq, h = 0.2, verbose = TRUE)
#' phyloseq::ntaxa(pq_glom)
#'
#' # Inspect the taxa-to-cluster mapping without merging
#' map <- phylo_glom_pq(pq, h = 0.2, return_map = TRUE)
#' head(map)
#' }
#'
#' \dontrun{
#' # Keep the taxonomy strictly conservative across merged taxa
#' pq_glom <- phylo_glom_pq(pq, h = 0.2, tax_adjust = 2)
#' }
#' @importFrom MiscMetabar verify_pq merge_taxa_vec
#' @export
phylo_glom_pq <- function(
  physeq,
  h,
  hclust_method = "average",
  tax_adjust = 1L,
  rank_propagation = TRUE,
  return_map = FALSE,
  verbose = FALSE
) {
  verify_pq(physeq)
  clusters <- phylo_clusters(physeq, h = h, hclust_method = hclust_method)

  if (return_map) {
    return(data.frame(
      taxa = names(clusters),
      cluster = unname(clusters),
      stringsAsFactors = FALSE
    ))
  }

  new_obj <- merge_taxa_vec(
    physeq,
    clusters,
    tax_adjust = tax_adjust,
    rank_propagation = rank_propagation
  )

  if (verbose) {
    cli::cli_inform(c(
      "v" = "{phyloseq::ntaxa(physeq)} taxa agglomerated into
             {phyloseq::ntaxa(new_obj)} taxa at h = {h}."
    ))
  }

  return(new_obj)
}

#' Scan several cophenetic thresholds before agglomerating
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Compute the number of taxa that [phylo_glom_pq()] would return
#' for a range of cophenetic thresholds, without building the agglomerated
#' objects. Useful to choose an informed value of `h` before running
#' [phylo_glom_pq()].
#'
#' @inheritParams phylo_glom_pq
#' @param h_values (required) Numeric vector of cophenetic distance thresholds
#'   to scan.
#'
#' @return A data.frame with one row per value of `h_values` and the columns
#'   `h` (the threshold), `n_taxa` (the resulting number of taxa) and
#'   `prop_taxa` (that number divided by the initial number of taxa).
#'
#' @author Adrien Taudière
#'
#' @seealso [phylo_glom_pq()]
#'
#' @examples
#' \donttest{
#' library(MiscMetabar)
#' data(data_fungi_mini)
#'
#' pq <- add_tree_pq(
#'   data_fungi_mini,
#'   use_taxo_to_build_tree = TRUE,
#'   compute_brlen = TRUE
#' )
#' phylo_glom_scan_pq(pq, h_values = seq(0, 1, by = 0.1))
#' }
#'
#' \dontrun{
#' scan <- phylo_glom_scan_pq(pq, h_values = seq(0, 1, by = 0.01))
#' plot(scan$h, scan$n_taxa, type = "l")
#' }
#' @importFrom MiscMetabar verify_pq
#' @export
phylo_glom_scan_pq <- function(
  physeq,
  h_values,
  hclust_method = "average"
) {
  verify_pq(physeq)
  n_taxa <- vapply(
    h_values,
    function(h) {
      length(unique(phylo_clusters(
        physeq,
        h = h,
        hclust_method = hclust_method
      )))
    },
    numeric(1)
  )

  data.frame(
    h = h_values,
    n_taxa = n_taxa,
    prop_taxa = n_taxa / phyloseq::ntaxa(physeq)
  )
}

#' Cut the cophenetic dendrogram of a phyloseq tree
#'
#' @description Internally used by [phylo_glom_pq()] and
#'   [phylo_glom_scan_pq()] to map taxa onto clusters.
#'
#' @inheritParams phylo_glom_pq
#' @return An integer vector of cluster memberships, named and ordered as
#'   `taxa_names(physeq)`.
#' @noRd
#' @keywords internal
phylo_clusters <- function(physeq, h, hclust_method = "average") {
  if (is.null(physeq@phy_tree)) {
    cli::cli_abort(c(
      "The {.arg physeq} object has no {.field phy_tree} slot.",
      "i" = "Use {.fn add_tree_pq} to attach a tree first."
    ))
  }

  tree <- phyloseq::phy_tree(physeq)
  if (!has_brlen(tree)) {
    cli::cli_abort(c(
      "The tree of {.arg physeq} has no branch length.",
      "i" = "Cophenetic distances are undefined without branch lengths.",
      "i" = "Use {.code add_tree_pq(physeq, compute_brlen = TRUE)} to add
             arbitrary ones."
    ))
  }

  if (!is.numeric(h) || length(h) != 1 || is.na(h) || h < 0) {
    cli::cli_abort("{.arg h} must be a single non-negative number.")
  }

  d <- stats::cophenetic(tree)
  taxa <- phyloseq::taxa_names(physeq)
  d <- d[taxa, taxa]

  hc <- stats::hclust(stats::as.dist(d), method = hclust_method)
  # Ties in the cophenetic matrix yield sub-epsilon height inversions
  hc$height <- cummax(hc$height)
  clusters <- stats::cutree(hc, h = h)

  clusters[taxa]
}
