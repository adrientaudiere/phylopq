#' Agglomerate taxa into network sequence clusters (NSC)
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Group the taxa of a \code{\link[phyloseq]{phyloseq-class}} object into
#' *network sequence clusters* (NSC), following the sequence-similarity network
#' (SSN) approach of Forster et al. (2019). All pairs of reference sequences are
#' compared with `vsearch --allpairs_global`, every pair at least `id` similar
#' becomes an edge of an undirected graph, and each connected component of that
#' graph is one cluster. Taxa of the same component are then merged with
#' [MiscMetabar::merge_taxa_vec()].
#'
#' A connected component is a *single-linkage* group: two sequences end up in
#' the same NSC as soon as a chain of pairwise similarities links them, even
#' when they are themselves below `id`. This is the point of the method — it
#' recovers the loose, chained variation of a species that a fixed centroid
#' threshold splits apart — but it also means an NSC can be much wider than
#' `id`, and that a single bridging sequence can fuse two clusters. Use
#' [ssn_glom_scan_pq()] to see how the number of clusters responds to `id`
#' before committing to a value.
#'
#' Unlike [MiscMetabar::postcluster_pq()], which assigns each sequence to a
#' centroid, and unlike [phylo_glom_pq()], which cuts a tree, `ssn_glom_pq()`
#' needs no centroid and no tree: only the `refseq` slot.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   with a non-empty `refseq` slot.
#' @param id Numeric in \[0, 1\], the pairwise-similarity threshold above which
#'   two sequences are linked by an edge. Default to 0.97.
#' @param iddef Integer, the identity definition used by vsearch: 0 (CD-HIT
#'   definition), 1 (edit distance, the default here and the definition used by
#'   Forster et al. 2019), 2 (edit distance excluding terminal gaps), 3
#'   (marine-biological definition) or 4 (edit distance excluding gaps of any
#'   kind).
#' @param vsearchpath Path to the vsearch executable. Default to NULL, in which
#'   case it is looked up with [MiscMetabar::find_vsearch()].
#' @param vsearch_args A length-one character vector of further arguments
#'   passed on to vsearch, e.g. `"--threads 4"`. Note that `--strand` is not a
#'   valid option of `--allpairs_global`. Default to `""`.
#' @param tax_adjust Handling of taxonomic disagreements within a cluster. See
#'   [MiscMetabar::merge_taxa_vec()]. Default to 1 (phyloseq-compatible).
#' @param rank_propagation Logical, default TRUE, whether bad ranks are
#'   propagated to lower ranks. See [MiscMetabar::merge_taxa_vec()].
#' @param return_map Logical, if TRUE return a data.frame mapping each taxon to
#'   its cluster instead of the agglomerated phyloseq object. Default to FALSE.
#' @param return_graph Logical, if TRUE return the `igraph` sequence-similarity
#'   network itself, with the cluster membership stored as the `nsc` vertex
#'   attribute. Default to FALSE. Takes precedence over `return_map`.
#' @param keep_temporary_files Logical, if TRUE the FASTA and `.uc` files
#'   handed to and produced by vsearch are kept. Default to FALSE.
#' @param verbose Logical, if TRUE report the number of edges, of clusters and
#'   of taxa before and after agglomeration. Default to FALSE.
#'
#' @return A \code{\link[phyloseq]{phyloseq-class}} object whose taxa are the
#'   network sequence clusters; or a data.frame with columns `taxa` and
#'   `cluster` when `return_map = TRUE`; or an `igraph` object when
#'   `return_graph = TRUE`.
#'
#' @details vsearch is an external program and is not installed by phylopq. Use
#'   [MiscMetabar::is_vsearch_installed()] to check for it and
#'   [MiscMetabar::install_vsearch()] to install it.
#'
#'   `--allpairs_global` aligns every pair of sequences, so the run is
#'   quadratic in the number of taxa. It is comfortable for a few thousand
#'   taxa and slow well beyond that; subset or pre-cluster first in that case.
#'
#'   Taxa that vsearch reports no hit for form their own single-taxon cluster
#'   rather than being dropped, so no abundance is lost by the agglomeration.
#'
#' @author Adrien Taudière
#'
#' @references
#' Forster D., Filker S., Kochems R., Breiner H.-W., Cordier T., Pawlowski J.,
#' Stoeck T. (2019) A comparison of different ciliate metabarcode genes as
#' bioindicators for environmental impact assessments of salmon aquaculture.
#' *Environmental Microbiology* 21(11):4429-4443.
#' \doi{10.1111/1462-2920.14764}
#'
#' @seealso [ssn_glom_scan_pq()], [phylo_glom_pq()],
#'   [MiscMetabar::postcluster_pq()], [MiscMetabar::merge_taxa_vec()]
#'
#' @examplesIf MiscMetabar::is_vsearch_installed() && requireNamespace("igraph", quietly = TRUE)
#' library(MiscMetabar)
#' data(data_fungi_mini)
#' df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 5000)
#' phyloseq::ntaxa(df)
#'
#' pq_nsc <- ssn_glom_pq(df, id = 0.9, verbose = TRUE)
#' phyloseq::ntaxa(pq_nsc)
#'
#' # Inspect the taxa-to-cluster mapping without merging
#' map <- ssn_glom_pq(df, id = 0.9, return_map = TRUE)
#' head(map)
#'
#' \dontrun{
#' # The network itself, to plot or to analyse further
#' net <- ssn_glom_pq(df, id = 0.9, return_graph = TRUE)
#' igraph::V(net)$nsc
#' plot(net, vertex.label = NA, vertex.size = 5)
#'
#' # Keep the taxonomy strictly conservative across merged taxa
#' pq_nsc <- ssn_glom_pq(df, id = 0.9, tax_adjust = 2)
#'
#' # An explicit path to the executable
#' pq_nsc <- ssn_glom_pq(df, id = 0.9, vsearchpath = "/usr/local/bin/vsearch")
#' }
#' @importFrom MiscMetabar verify_pq merge_taxa_vec
#' @export
ssn_glom_pq <- function(
  physeq,
  id = 0.97,
  iddef = 1,
  vsearchpath = NULL,
  vsearch_args = "",
  tax_adjust = 1L,
  rank_propagation = TRUE,
  return_map = FALSE,
  return_graph = FALSE,
  keep_temporary_files = FALSE,
  verbose = FALSE
) {
  verify_pq(physeq)
  validate_ssn_id(id)

  pairs <- ssn_pairs(
    physeq,
    id = id,
    iddef = iddef,
    vsearchpath = vsearchpath,
    vsearch_args = vsearch_args,
    keep_temporary_files = keep_temporary_files,
    verbose = verbose
  )

  taxa <- phyloseq::taxa_names(physeq)
  clusters <- ssn_clusters(pairs, taxa = taxa, id = id)

  if (verbose) {
    cli::cli_inform(c(
      "i" = "{nrow(pairs[pairs$id >= id, , drop = FALSE])} edge{?s} at
             id = {id} link{?s/} {length(taxa)} taxa into
             {length(unique(clusters))} cluster{?s}."
    ))
  }

  if (return_graph) {
    return(ssn_graph(pairs, taxa = taxa, id = id, clusters = clusters))
  }

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
             {phyloseq::ntaxa(new_obj)} network sequence cluster{?s} at
             id = {id}."
    ))
  }

  return(new_obj)
}

#' Scan several similarity thresholds before agglomerating
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Compute the number of network sequence clusters that [ssn_glom_pq()] would
#' return for a range of similarity thresholds, without building the
#' agglomerated objects. Useful to choose an informed value of `id`, as a
#' single-linkage network is sensitive to it.
#'
#' vsearch is run **once**, at the lowest threshold of `id_values`, and the
#' resulting pairwise identities are then filtered at each threshold in turn —
#' so scanning twenty thresholds costs one alignment run, not twenty.
#'
#' @inheritParams ssn_glom_pq
#' @param id_values (required) Numeric vector of similarity thresholds in
#'   \[0, 1\] to scan.
#'
#' @return A data.frame with one row per value of `id_values` and the columns
#'   `id` (the threshold), `n_edges` (the number of pairs at or above it),
#'   `n_taxa` (the resulting number of clusters) and `prop_taxa` (that number
#'   divided by the initial number of taxa).
#'
#' @author Adrien Taudière
#'
#' @seealso [ssn_glom_pq()], [phylo_glom_scan_pq()]
#'
#' @examplesIf MiscMetabar::is_vsearch_installed() && requireNamespace("igraph", quietly = TRUE)
#' library(MiscMetabar)
#' data(data_fungi_mini)
#' df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 5000)
#'
#' ssn_glom_scan_pq(df, id_values = seq(0.8, 1, by = 0.05))
#'
#' \dontrun{
#' scan <- ssn_glom_scan_pq(df, id_values = seq(0.7, 1, by = 0.01))
#' plot(scan$id, scan$n_taxa, type = "l", xlab = "id", ylab = "Nb of NSC")
#' }
#' @importFrom MiscMetabar verify_pq
#' @export
ssn_glom_scan_pq <- function(
  physeq,
  id_values,
  iddef = 1,
  vsearchpath = NULL,
  vsearch_args = "",
  keep_temporary_files = FALSE,
  verbose = FALSE
) {
  verify_pq(physeq)
  if (length(id_values) == 0) {
    cli::cli_abort("{.arg id_values} must hold at least one threshold.")
  }
  lapply(id_values, validate_ssn_id)

  pairs <- ssn_pairs(
    physeq,
    id = min(id_values),
    iddef = iddef,
    vsearchpath = vsearchpath,
    vsearch_args = vsearch_args,
    keep_temporary_files = keep_temporary_files,
    verbose = verbose
  )

  taxa <- phyloseq::taxa_names(physeq)
  n_edges <- vapply(
    id_values,
    function(id) {
      sum(pairs$id >= id)
    },
    numeric(1)
  )
  n_taxa <- vapply(
    id_values,
    function(id) {
      length(unique(ssn_clusters(pairs, taxa = taxa, id = id)))
    },
    numeric(1)
  )

  data.frame(
    id = id_values,
    n_edges = n_edges,
    n_taxa = n_taxa,
    prop_taxa = n_taxa / length(taxa)
  )
}

#' Check the `id` argument of [ssn_glom_pq()]
#'
#' @inheritParams ssn_glom_pq
#' @return `id`, invisibly, or an error.
#' @noRd
#' @keywords internal
validate_ssn_id <- function(id) {
  if (!is.numeric(id) || length(id) != 1 || is.na(id) || id < 0 || id > 1) {
    cli::cli_abort(
      "{.arg id} must be a single number between 0 and 1, not {.val {id}}."
    )
  }
  invisible(id)
}

#' All pairs of reference sequences at least `id` similar
#'
#' Runs `vsearch --allpairs_global` on the `refseq` slot and reads back the
#' `.uc` file it writes. Only `H` records carry an identity; the `N` records of
#' sequences without any hit are dropped here and recovered as single-taxon
#' clusters by [ssn_clusters()].
#'
#' @inheritParams ssn_glom_pq
#' @return A data.frame with columns `from`, `to` (taxa names) and `id` (the
#'   pairwise identity as a fraction).
#' @noRd
#' @keywords internal
ssn_pairs <- function(
  physeq,
  id,
  iddef = 1,
  vsearchpath = NULL,
  vsearch_args = "",
  keep_temporary_files = FALSE,
  verbose = FALSE
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg igraph} is required by {.fn ssn_glom_pq}.",
      "i" = "Install it with {.code install.packages(\"igraph\")}."
    ))
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg Biostrings} is required by {.fn ssn_glom_pq}."
    )
  }
  if (is.null(physeq@refseq)) {
    cli::cli_abort(c(
      "The {.arg physeq} object has no {.field refseq} slot.",
      "i" = "{.fn ssn_glom_pq} builds its network from reference sequences."
    ))
  }
  if (phyloseq::ntaxa(physeq) < 2) {
    cli::cli_abort(
      "At least two taxa are required to build a sequence-similarity network,
       not {phyloseq::ntaxa(physeq)}."
    )
  }

  exec <- resolve_vsearch_exec(vsearchpath)

  fasta_file <- tempfile("phylopq_ssn_", fileext = ".fasta")
  uc_file <- tempfile("phylopq_ssn_", fileext = ".uc")
  if (!keep_temporary_files) {
    on.exit(unlink(c(fasta_file, uc_file)), add = TRUE)
  }
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(physeq@refseq),
    fasta_file
  )

  # `--notrunclabels` keeps taxa names holding whitespace intact, so that the
  # labels read back from the `.uc` file still match `taxa_names(physeq)`.
  log_file <- tempfile("phylopq_ssn_", fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)
  status <- system2(
    exec,
    paste(
      "--allpairs_global",
      shQuote(fasta_file),
      "--uc",
      shQuote(uc_file),
      "--iddef",
      iddef,
      "--id",
      id,
      "--notrunclabels",
      vsearch_args,
      if (verbose) "" else "--quiet"
    ),
    stdout = if (verbose) "" else FALSE,
    stderr = log_file
  )

  if (status != 0 || !file.exists(uc_file)) {
    msg <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else ""
    msg <- utils::head(msg[nzchar(msg)], 5)
    cli::cli_abort(c(
      "{.field vsearch} exited with status {status} and wrote no
       {.field .uc} file.",
      if (length(msg) > 0) {
        stats::setNames(msg, rep("x", length(msg)))
      },
      "i" = "Check that {.path {exec}} runs from a terminal.",
      "i" = "Use {.code verbose = TRUE} to see its output."
    ))
  }

  read_uc_pairs(uc_file)
}

#' Read the `H` records of a vsearch `.uc` file
#'
#' The `.uc` format is tab-separated and column-positional: column 1 is the
#' record type, column 4 the identity as a percentage, and columns 9 and 10 the
#' query and target labels. Non-hit (`N`) records carry `*` in place of the
#' identity and of the target.
#'
#' @param path Path to a `.uc` file.
#' @return A data.frame with columns `from`, `to` and `id` (a fraction).
#' @noRd
#' @keywords internal
read_uc_pairs <- function(path) {
  empty <- data.frame(
    from = character(0),
    to = character(0),
    id = numeric(0),
    stringsAsFactors = FALSE
  )

  if (file.size(path) == 0) {
    return(empty)
  }

  uc <- utils::read.table(
    path,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    na.strings = "*"
  )

  hits <- uc[uc[[1]] == "H" & !is.na(uc[[4]]), , drop = FALSE]
  if (nrow(hits) == 0) {
    return(empty)
  }

  data.frame(
    from = as.character(hits[[9]]),
    to = as.character(hits[[10]]),
    id = as.numeric(hits[[4]]) / 100,
    stringsAsFactors = FALSE
  )
}

#' Connected components of the sequence-similarity network
#'
#' @param pairs A data.frame as returned by [ssn_pairs()].
#' @param taxa Character vector of every taxa name, so that taxa without any
#'   edge become their own single-taxon cluster instead of vanishing.
#' @inheritParams ssn_glom_pq
#' @return An integer vector of cluster memberships, named and ordered as
#'   `taxa`.
#' @noRd
#' @keywords internal
ssn_clusters <- function(pairs, taxa, id) {
  net <- ssn_graph(pairs, taxa = taxa, id = id)
  clusters <- igraph::components(net)$membership
  clusters[taxa]
}

#' Build the sequence-similarity network at one threshold
#'
#' @inheritParams ssn_clusters
#' @param clusters Optional membership vector to store as the `nsc` vertex
#'   attribute, as returned by [ssn_clusters()].
#' @return An undirected `igraph` object with one vertex per taxa name and one
#'   edge per pair at or above `id`, carrying the identity as the `id` edge
#'   attribute.
#' @noRd
#' @keywords internal
ssn_graph <- function(pairs, taxa, id, clusters = NULL) {
  kept <- pairs[pairs$id >= id, , drop = FALSE]

  unknown <- setdiff(c(kept$from, kept$to), taxa)
  if (length(unknown) > 0) {
    cli::cli_abort(c(
      "{.field vsearch} returned {length(unknown)} label{?s} that
       {?is/are} not {?a taxa name/taxa names}: {.val {utils::head(unknown)}}.",
      "i" = "Taxa names holding whitespace are truncated in a FASTA header.",
      "i" = "Rename the taxa, e.g. with
             {.code taxa_names(physeq) <- paste0(\"ASV_\", seq_len(ntaxa(physeq)))}."
    ))
  }

  net <- igraph::graph_from_data_frame(
    kept,
    directed = FALSE,
    vertices = data.frame(name = taxa, stringsAsFactors = FALSE)
  )

  if (!is.null(clusters)) {
    net <- igraph::set_vertex_attr(
      net,
      "nsc",
      value = unname(clusters[igraph::V(net)$name])
    )
  }

  net
}

#' Resolve the path to the vsearch executable
#'
#' @inheritParams ssn_glom_pq
#' @return A length-one character path to the executable.
#' @noRd
#' @keywords internal
resolve_vsearch_exec <- function(vsearchpath) {
  if (is.null(vsearchpath)) {
    vsearchpath <- MiscMetabar::find_vsearch()
  }
  if (!MiscMetabar::is_vsearch_installed(vsearchpath)) {
    cli::cli_abort(c(
      "The {.field vsearch} executable was not found.",
      "i" = "Install it with {.code MiscMetabar::install_vsearch()}, or e.g.
             {.code sudo apt install vsearch} or
             {.code conda install -c bioconda vsearch}.",
      "i" = "Give its path with {.arg vsearchpath}, or set it once with
             {.code options(MiscMetabar.vsearchpath = \"/path/to/vsearch\")}.",
      "i" = "Installation instructions:
             {.url https://github.com/torognes/vsearch}."
    ))
  }
  vsearchpath
}
