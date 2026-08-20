#' Compare several species delimitations of a phyloseq object
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Run ABGD and ASAP several times on the reference sequences of a
#' \code{\link[phyloseq]{phyloseq-class}} object, join the resulting partitions
#' into a single table with [delimtools::delim_join()], and display them side
#' by side along a phylogenetic tree with [delimtools::delim_autoplot()].
#'
#' Where [delim_pq()] answers "which species does this taxon belong to?",
#' `delim_multi_pq()` answers "how stable is that answer?". A single barcode-gap
#' run is very sensitive to its parameters, so the usual practice is to scan a
#' range of relative gap widths and to compare the methods with each other
#' before trusting any one partition.
#'
#' `slopes` drives the ABGD runs — one run per value, since the relative gap
#' width is the parameter ABGD is most sensitive to. ASAP has no equivalent
#' parameter (it scores and ranks partitions itself), so it contributes a
#' single run whatever `slopes` contains.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   with a non-empty `refseq` slot.
#' @param methods Character vector of delimitation programs to run, any of
#'   `"abgd"` and `"asap"`. Default to both.
#' @param slopes Numeric vector of relative gap widths, one ABGD run per value.
#'   Default to `c(0.5, 1, 1.5)`. Ignored when `"abgd"` is not in `methods`.
#' @param model Integer, the evolutionary model used by the external programs.
#'   0: Kimura-2P, 1: Jukes-Cantor, 2: Tamura-Nei, 3: simple p-distance
#'   (default).
#' @param exe Named character vector of paths to the executables, e.g.
#'   `c(abgd = "/usr/local/bin/abgd")`. Default to NULL, in which case each
#'   program is looked up with [is_delim_installed()].
#' @param tree The tree used to order and display the partitions. Either a
#'   `phylo` object (`ape` package), a `treedata` object (`tidytree` package)
#'   or NULL (default). When NULL the `phy_tree` slot of `physeq` is used, and
#'   failing that a tree is built with [MiscMetabar::build_phytree_pq()] if
#'   `build_tree = TRUE`.
#' @param build_tree Logical, if TRUE and no tree is available a maximum
#'   likelihood tree is built from the alignment with
#'   [MiscMetabar::build_phytree_pq()]. Default to FALSE, because building a
#'   tree is by far the slowest step of this function.
#' @param plot Logical, if TRUE (default) a [delimtools::delim_autoplot()]
#'   figure is returned in the `plot` element. Requires a tree.
#' @param consensus Logical, if TRUE (default) a consensus column computed by
#'   [delimtools::delim_consensus()] is added to the figure. Silently turned
#'   off when fewer than two delimitations succeeded.
#' @param n_match Integer, the number of delimitations that must agree for the
#'   consensus to call a species. Default to NULL, i.e. a simple majority
#'   (at least two).
#' @param delim_order Character vector giving the order of the columns in the
#'   figure. Default to NULL, i.e. the order the runs were performed in.
#' @param tbl_labs A two-column data.frame mapping tip labels (`label`) to the
#'   labels shown on the figure (`labs`). Default to NULL, i.e. the tip labels
#'   themselves.
#' @param col_vec Character vector of colours passed to
#'   [delimtools::delim_autoplot()]. Default to NULL, i.e. a palette built with
#'   [delimtools::delim_brewer()].
#' @param widths Relative widths of the tree panel and of the partition panel,
#'   passed to [delimtools::delim_autoplot()]. Default to `c(0.5, 0.3)`, a
#'   little wider than the `delimtools` default so that several runs fit.
#' @param label_angle Angle of the run names above the partition panel.
#'   Default to 45 degrees, because run names such as `abgd_slope0.5` overlap
#'   when written horizontally. Set to 0 for the `delimtools` look.
#' @param align Logical, if TRUE (default) reference sequences of unequal
#'   length are aligned with [align_pq()] before being submitted. The alignment
#'   is computed once and reused by every run.
#' @param align_method Aligner used when `align = TRUE`, either `"decipher"`
#'   (default, pure R) or `"mafft"` (external program, much faster on large
#'   `refseq` slots). See [align_pq()].
#' @param mafft_exec Path to the MAFFT executable. Only used when
#'   `align_method = "mafft"`. Default to NULL, i.e. the usual lookup of
#'   [is_mafft_installed()].
#' @param outfolder Path to the folder the external programs write into. Each
#'   run gets its own sub-folder, so that a run producing nothing cannot
#'   inherit another run's partitions. Default to NULL, i.e. a temporary
#'   location.
#' @param keep_temporary_files Logical, if TRUE the FASTA file and the output
#'   folders are kept. Default to FALSE.
#' @param verbose Logical, if TRUE report each run and print the
#'   [delimtools::report_delim()] summary. Default to FALSE.
#' @param ... Additional arguments passed on to
#'   [MiscMetabar::build_phytree_pq()] when `build_tree = TRUE`.
#'
#' @return A list with five elements:
#' \describe{
#'   \item{delim}{A named list of the partition tables, one per run.}
#'   \item{joined}{A tibble with a `labels` column and one column per run,
#'     whose partition names are harmonised across runs by
#'     [delimtools::delim_join()], so that the same partition carries the same
#'     name everywhere.}
#'   \item{summary}{A data.frame giving the number of partitions each run
#'     found.}
#'   \item{plot}{The [delimtools::delim_autoplot()] figure, or NULL when
#'     `plot = FALSE` or no tree was available.}
#'   \item{tree}{The `phylo` object used for the figure, or NULL.}
#' }
#'
#' @details Run names are built from the method and its parameter, e.g.
#'   `abgd_slope0.5`. Internally each run is submitted to
#'   [delimtools::delim_join()] under a digit-free alias, because that function
#'   strips every digit from the delimitation names before matching partitions
#'   across runs, which would collapse `abgd_slope0.5` and `abgd_slope1.5` into
#'   a single column. The user-facing names are restored afterwards.
#'
#'   ABGD and ASAP are not available on Windows, and both need an aligned
#'   FASTA — hence the `align = TRUE` default.
#'
#' @author Adrien Taudière
#'
#' @references
#' Puillandre N., Lambert A., Brouillet S., Achaz G. (2012) ABGD, Automatic
#' Barcode Gap Discovery for primary species delimitation. *Molecular Ecology*
#' 21(8):1864-1877. \doi{10.1111/j.1365-294X.2011.05239.x}
#'
#' Puillandre N., Brouillet S., Achaz G. (2021) ASAP: assemble species by
#' automatic partitioning. *Molecular Ecology Resources* 21:609-620.
#' \doi{10.1111/1755-0998.13281}
#'
#' @seealso [delim_pq()], [is_delim_installed()], [align_pq()],
#'   [MiscMetabar::build_phytree_pq()], [delimtools::delim_autoplot()]
#'
#' @examplesIf phylopq::is_delim_installed("asap") && phylopq::is_delim_installed("abgd")
#' \donttest{
#' library(MiscMetabar)
#' data(data_fungi_mini)
#' df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 9000)
#'
#' # Scan three ABGD gap widths, add one ASAP run, and order the result along
#' # a taxonomy-based tree
#' res <- delim_multi_pq(
#'   df,
#'   slopes = c(0.5, 1, 1.5),
#'   tree = taxo2tree(df),
#'   verbose = TRUE
#' )
#' res$summary
#' res$joined
#' res$plot
#' }
#'
#' \dontrun{
#' # Without a tree, only the tables are returned
#' res <- delim_multi_pq(df, plot = FALSE)
#'
#' # Build a maximum likelihood tree on the fly (slow) and align with MAFFT
#' res <- delim_multi_pq(
#'   df,
#'   build_tree = TRUE,
#'   nb_bootstrap = 5,
#'   align_method = "mafft"
#' )
#'
#' # ABGD only, on a finer grid of gap widths
#' res <- delim_multi_pq(df, methods = "abgd", slopes = seq(0.4, 1.6, by = 0.2))
#' }
#' @importFrom MiscMetabar verify_pq
#' @export
delim_multi_pq <- function(
  physeq,
  methods = c("abgd", "asap"),
  slopes = c(0.5, 1, 1.5),
  model = 3,
  exe = NULL,
  tree = NULL,
  build_tree = FALSE,
  plot = TRUE,
  consensus = TRUE,
  n_match = NULL,
  delim_order = NULL,
  tbl_labs = NULL,
  col_vec = NULL,
  widths = c(0.5, 0.3),
  label_angle = 45,
  align = TRUE,
  align_method = c("decipher", "mafft"),
  mafft_exec = NULL,
  outfolder = NULL,
  keep_temporary_files = FALSE,
  verbose = FALSE,
  ...
) {
  verify_pq(physeq)
  methods <- match.arg(methods, c("abgd", "asap"), several.ok = TRUE)
  align_method <- match.arg(align_method)

  if (!requireNamespace("delimtools", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg delimtools} is required by {.fn delim_multi_pq}."
    )
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg Biostrings} is required by {.fn delim_multi_pq}."
    )
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg tibble} is required by {.fn delim_multi_pq}."
    )
  }
  if (is.null(physeq@refseq)) {
    cli::cli_abort(c(
      "The {.arg physeq} object has no {.field refseq} slot.",
      "i" = "{.fn delim_multi_pq} delimits species from reference sequences."
    ))
  }

  runs <- delim_run_grid(methods, slopes)

  # Align once and share the FASTA: the alignment is the same for every run
  # and is often the most expensive step after tree building.
  dna <- align_refseq(
    Biostrings::DNAStringSet(physeq@refseq),
    align,
    align_method,
    mafft_exec
  )
  fasta_file <- file.path(tempdir(), "delim_multi_pq.fasta")
  Biostrings::writeXStringSet(dna, fasta_file)
  if (!keep_temporary_files) {
    on.exit(unlink(fasta_file), add = TRUE)
  }

  if (is.null(outfolder)) {
    outfolder <- tempfile("phylopq_delim_multi_")
    if (!keep_temporary_files) {
      on.exit(unlink(outfolder, recursive = TRUE), add = TRUE)
    }
  }
  dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)

  exes <- resolve_delim_exes(methods, exe)

  delim <- stats::setNames(vector("list", nrow(runs)), runs$label)
  for (i in seq_len(nrow(runs))) {
    delim[[i]] <- delim_one_run(
      run = runs[i, ],
      fasta_file = fasta_file,
      exes = exes,
      model = model,
      outfolder = outfolder,
      verbose = verbose
    )
  }

  failed <- vapply(delim, is.null, logical(1))
  if (all(failed)) {
    cli::cli_abort(c(
      "None of the {nrow(runs)} run{?s} returned a partition.",
      "i" = "See the {.pkg delimtools} messages above for the reason."
    ))
  }
  if (any(failed)) {
    cli::cli_warn(
      "{sum(failed)} run{?s} returned no partition and {?is/are} dropped:
       {.val {names(delim)[failed]}}."
    )
    delim <- delim[!failed]
    runs <- runs[!failed, , drop = FALSE]
  }

  joined <- delim_join_runs(delim, runs)

  summary_df <- data.frame(
    run = setdiff(colnames(joined), "labels"),
    n_partitions = vapply(
      setdiff(colnames(joined), "labels"),
      function(nm) {
        length(unique(joined[[nm]]))
      },
      integer(1)
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    invisible(delimtools::report_delim(joined))
  }

  tree <- resolve_delim_tree(physeq, tree, build_tree, dna, plot, verbose, ...)

  p <- NULL
  if (plot && !is.null(tree)) {
    p <- delim_autoplot_pq(
      joined = joined,
      tree = tree,
      consensus = consensus,
      n_match = n_match,
      delim_order = delim_order,
      tbl_labs = tbl_labs,
      col_vec = col_vec,
      widths = widths,
      label_angle = label_angle
    )
  }

  return(list(
    delim = delim,
    joined = joined,
    summary = summary_df,
    plot = p,
    tree = tree
  ))
}

#' Build the grid of delimitation runs
#'
#' ABGD gets one run per slope; ASAP has no relative gap width, so it gets a
#' single run whatever `slopes` contains.
#'
#' Aliases must be free of digits: [delimtools::delim_join()] strips every
#' digit from the delimitation name before matching partitions across runs, so
#' `abgd_slope0.5` and `abgd_slope1.5` would both become `abgd_slope.` and
#' collide.
#'
#' @inheritParams delim_multi_pq
#' @return A data.frame with columns `method`, `slope`, `label` and `alias`.
#' @noRd
#' @keywords internal
delim_run_grid <- function(methods, slopes) {
  if ("abgd" %in% methods) {
    if (length(slopes) == 0 || !is.numeric(slopes) || anyNA(slopes)) {
      cli::cli_abort(
        "{.arg slopes} must be a non-empty numeric vector without missing
         values."
      )
    }
    slopes <- unique(slopes)
  }

  grid <- do.call(
    rbind,
    lapply(methods, function(m) {
      if (identical(m, "abgd")) {
        data.frame(
          method = "abgd",
          slope = slopes,
          # Format value by value: `format()` on the whole vector pads to a
          # common width, so the name of a run would depend on the other
          # slopes requested in the same call.
          label = paste0(
            "abgd_slope",
            vapply(slopes, format, character(1), trim = TRUE)
          ),
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          method = "asap",
          slope = NA_real_,
          label = "asap",
          stringsAsFactors = FALSE
        )
      }
    })
  )

  grid$alias <- paste0(grid$method, alias_letters(nrow(grid)))
  grid
}

#' Digit-free aliases, in the spirit of Excel column names
#'
#' @param n Number of aliases to generate.
#' @return A character vector of length `n`: A, B, ..., Z, AA, AB, ...
#' @noRd
#' @keywords internal
alias_letters <- function(n) {
  out <- character(n)
  for (i in seq_len(n)) {
    idx <- i
    label <- ""
    while (idx > 0) {
      rest <- (idx - 1) %% 26
      label <- paste0(LETTERS[rest + 1], label)
      idx <- (idx - 1) %/% 26
    }
    out[i] <- label
  }
  out
}

#' Resolve the executables needed by [delim_multi_pq()]
#'
#' @inheritParams delim_multi_pq
#' @return A named character vector of paths, one per method.
#' @noRd
#' @keywords internal
resolve_delim_exes <- function(methods, exe) {
  stats::setNames(
    vapply(
      methods,
      function(m) {
        resolve_delim_exe(m, if (is.null(exe)) NULL else unname(exe[m]))
      },
      character(1)
    ),
    methods
  )
}

#' Perform one delimitation run of [delim_multi_pq()]
#'
#' @inheritParams delim_multi_pq
#' @param run A one-row data.frame of the run grid.
#' @param fasta_file Path to the shared aligned FASTA file.
#' @param exes A named character vector of executable paths.
#' @return A `tbl_df` with columns `labels` and the run alias, or NULL when the
#'   program returned nothing.
#' @noRd
#' @keywords internal
delim_one_run <- function(run, fasta_file, exes, model, outfolder, verbose) {
  run_folder <- file.path(outfolder, run$alias)
  dir.create(run_folder, recursive = TRUE, showWarnings = FALSE)

  if (verbose) {
    cli::cli_inform(c("i" = "Running {.field {run$label}}."))
  }

  res <- tryCatch(
    if (identical(run$method, "abgd")) {
      quiet_unless(
        delimtools::abgd_tbl(
          infile = fasta_file,
          exe = exes[["abgd"]],
          slope = run$slope,
          model = model,
          outfolder = run_folder,
          webserver = NULL,
          delimname = run$alias
        ),
        verbose
      )
    } else {
      asap_tbl_safely(
        infile = fasta_file,
        exe = exes[["asap"]],
        model = model,
        outfolder = run_folder,
        webserver = NULL,
        delimname = run$alias,
        verbose = verbose
      )
    },
    error = function(e) {
      cli::cli_warn(c(
        "{.field {run$label}} failed and is dropped.",
        "x" = conditionMessage(e)
      ))
      NULL
    }
  )

  res
}

#' Evaluate an expression, hiding its messages unless asked
#'
#' `delimtools::abgd_tbl()` reports the folder it wrote into on every call,
#' which is noise when a dozen runs are chained.
#'
#' @param expr An expression.
#' @param verbose Logical, if TRUE messages are let through.
#' @return The value of `expr`.
#' @noRd
#' @keywords internal
quiet_unless <- function(expr, verbose) {
  if (verbose) {
    return(expr)
  }
  suppressMessages(expr)
}

#' Join the runs of [delim_multi_pq()] under digit-free aliases
#'
#' @inheritParams delim_multi_pq
#' @param delim A named list of partition tables.
#' @param runs The run grid, restricted to the runs that succeeded.
#' @return A tibble with a `labels` column and one column per run, named after
#'   the user-facing run labels.
#' @noRd
#' @keywords internal
delim_join_runs <- function(delim, runs) {
  # `delim_join()` aborts on any NA, so keep only the labels every run
  # delimited rather than letting it fail on a partial run.
  shared <- Reduce(
    intersect,
    lapply(delim, function(x) {
      x$labels[!is.na(x[[2]])]
    })
  )
  if (length(shared) == 0) {
    cli::cli_abort(
      "No sequence was delimited by every run, so the partitions cannot be
       compared."
    )
  }
  dropped <- setdiff(delim[[1]]$labels, shared)
  if (length(dropped) > 0) {
    cli::cli_warn(
      "{length(dropped)} taxa {?was/were} not delimited by every run and
       {?is/are} dropped."
    )
  }
  delim <- lapply(delim, function(x) {
    x[match(shared, x$labels), , drop = FALSE]
  })

  joined <- if (length(delim) == 1) {
    as.data.frame(delim[[1]], stringsAsFactors = FALSE)
  } else {
    suppressMessages(delimtools::delim_join(unname(delim)))
  }
  joined <- tibble::as_tibble(joined)

  # A single run never goes through `delim_join()`, so its partition column is
  # still the raw integer the program returned; name it like the others.
  if (length(delim) == 1) {
    joined[[2]] <- paste0("sp", joined[[2]])
  }

  aliases <- colnames(joined)[-1]
  colnames(joined)[-1] <- runs$label[match(aliases, runs$alias)]
  joined
}

#' Resolve the tree used by [delim_multi_pq()]
#'
#' @inheritParams delim_multi_pq
#' @param dna The aligned `DNAStringSet`, used when a tree must be built.
#' @return A `phylo` object, or NULL when no tree is available.
#' @noRd
#' @keywords internal
resolve_delim_tree <- function(
  physeq,
  tree,
  build_tree,
  dna,
  plot,
  verbose,
  ...
) {
  if (inherits(tree, "treedata")) {
    return(tree@phylo)
  }
  if (inherits(tree, "phylo")) {
    return(tree)
  }
  if (!is.null(tree)) {
    cli::cli_abort(
      "{.arg tree} must be a {.cls phylo} or {.cls treedata} object, not
       {.cls {class(tree)}}."
    )
  }

  if (!is.null(physeq@phy_tree)) {
    if (verbose) {
      cli::cli_inform(c(
        "i" = "Using the {.field phy_tree} slot of
                              {.arg physeq}."
      ))
    }
    return(physeq@phy_tree)
  }

  if (build_tree) {
    if (verbose) {
      cli::cli_inform(c(
        "i" = "Building a tree with
                              {.fn MiscMetabar::build_phytree_pq}."
      ))
    }
    physeq_ali <- physeq
    physeq_ali@refseq <- dna
    return(MiscMetabar::build_phytree_pq(physeq_ali, ...)$ML$tree)
  }

  if (plot) {
    cli::cli_warn(c(
      "No tree available, so no figure is produced.",
      "i" = "Give one with {.arg tree}, add one to {.arg physeq} with
             {.fn add_tree_pq}, or use {.code build_tree = TRUE}.",
      "i" = "Use {.code plot = FALSE} to silence this warning."
    ))
  }
  NULL
}

#' Draw the [delimtools::delim_autoplot()] figure of [delim_multi_pq()]
#'
#' `delim_autoplot()` needs a `treedata` object carrying a `posterior` column
#' (ultrametric trees) or a `support` column (all others); a plain `phylo`
#' converted with [tidytree::as.treedata()] has neither, and the figure fails
#' while computing its aesthetics. Both columns are therefore filled from the
#' node labels when those are numeric, and with zeros otherwise, so that no
#' support point is drawn rather than the figure erroring.
#'
#' @inheritParams delim_multi_pq
#' @param joined The joined partition table.
#' @param tree A `phylo` object.
#' @return A `patchwork` object.
#' @noRd
#' @keywords internal
delim_autoplot_pq <- function(
  joined,
  tree,
  consensus,
  n_match,
  delim_order,
  tbl_labs,
  col_vec,
  widths,
  label_angle
) {
  pkgs <- c("ggplot2", "ggtree", "patchwork", "tidytree", "treeio", "tibble")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cli::cli_abort(c(
        "Package {.pkg {pkg}} is required to draw the delimitation figure.",
        "i" = "Use {.code plot = FALSE} to skip it."
      ))
    }
  }

  shared <- intersect(tree$tip.label, joined$labels)
  if (length(shared) == 0) {
    cli::cli_abort(c(
      "No tip label of the tree matches a delimited sequence.",
      "i" = "First tips: {.val {utils::head(tree$tip.label, 3)}}.",
      "i" = "First labels: {.val {utils::head(joined$labels, 3)}}."
    ))
  }
  if (length(shared) < length(joined$labels)) {
    cli::cli_warn(
      "{length(joined$labels) - length(shared)} delimited sequence{?s}
       {?is/are} absent from the tree and {?is/are} not shown."
    )
  }
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, shared))
  joined <- joined[joined$labels %in% shared, , drop = FALSE]

  if (length(tree$edge.length) == 0) {
    tree <- ape::compute.brlen(tree)
  }

  runs <- setdiff(colnames(joined), "labels")
  if (is.null(delim_order)) {
    delim_order <- runs
  }
  if (consensus && length(runs) < 2) {
    consensus <- FALSE
  }
  if (consensus && is.null(n_match)) {
    n_match <- max(2, ceiling(length(runs) / 2))
  }
  if (is.null(tbl_labs)) {
    tbl_labs <- tibble::tibble(label = tree$tip.label, labs = tree$tip.label)
  }
  if (is.null(col_vec)) {
    col_vec <- suppressWarnings(delimtools::delim_brewer(
      joined,
      package = NULL,
      palette = NULL,
      seed = NULL
    ))
  }

  p <- delimtools::delim_autoplot(
    delim = joined,
    tr = as_delim_treedata(tree),
    consensus = consensus,
    n_match = n_match,
    delim_order = delim_order,
    tbl_labs = tbl_labs,
    col_vec = col_vec,
    widths = widths
  )

  # The second panel of the patchwork holds the partition tiles; its run names
  # are written horizontally by `delim_autoplot()` and overlap as soon as a few
  # runs are compared.
  if (!is.null(label_angle) && label_angle != 0 && length(p) == 2) {
    p[[2]] <- p[[2]] +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = label_angle,
          hjust = 0,
          vjust = 0
        )
      )
  }

  p
}

#' Convert a `phylo` object to the `treedata` object [delimtools::delim_autoplot()] expects
#'
#' @param tree A `phylo` object.
#' @return A `treedata` object with `support` and `posterior` node data.
#' @noRd
#' @keywords internal
as_delim_treedata <- function(tree) {
  n_tip <- ape::Ntip(tree)
  node_support <- suppressWarnings(as.numeric(tree$node.label))
  if (length(node_support) != tree$Nnode) {
    node_support <- rep(NA_real_, tree$Nnode)
  }
  node_support[is.na(node_support)] <- 0
  support <- c(rep(0, n_tip), node_support)

  # `delim_autoplot()` reads bootstrap support on a 0-100 scale and posterior
  # probabilities on a 0-1 scale; node labels carry only one of the two.
  posterior <- ifelse(support > 1, support / 100, support)

  td <- treeio::as.treedata(tree)
  td@data <- tibble::tibble(
    node = seq_len(n_tip + tree$Nnode),
    support = support,
    posterior = posterior
  )
  td
}
