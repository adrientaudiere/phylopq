#' Species delimitation of a phyloseq object with ABGD or ASAP
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Run a single-locus species-delimitation analysis on the
#' reference sequences of a \code{\link[phyloseq]{phyloseq-class}} object and
#' map the resulting partitions back onto its taxa.
#'
#' Two barcode-gap methods are available through the `delimtools` package:
#' ABGD (Automatic Barcode Gap Discovery, Puillandre et al. 2012) and ASAP
#' (Assemble Species by Automatic Partitioning, Puillandre et al. 2021). Both
#' are external command-line programs: either give the path to the executable
#' with `exe`, or pass a result file downloaded from the corresponding web
#' server with `webserver`.
#'
#' By default taxa belonging to the same partition are merged with
#' [MiscMetabar::merge_taxa_vec()], so `delim_pq()` acts as a post-clustering
#' step comparable to [MiscMetabar::postcluster_pq()], but driven by a
#' barcode gap rather than by a fixed similarity threshold.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   with a non-empty `refseq` slot.
#' @param method One of `"asap"` (default) or `"abgd"`.
#' @param exe Path to the ABGD or ASAP executable. Ignored when `webserver` is
#'   given. Default to NULL, in which case the executable is looked up with
#'   [is_delim_installed()]: first the `phylopq.asappath` / `phylopq.abgdpath`
#'   option, then a copy installed by [install_asap()] / [install_abgd()],
#'   then the system `PATH`.
#' @param model Integer, the evolutionary model used by the external program.
#'   0: Kimura-2P, 1: Jukes-Cantor, 2: Tamura-Nei, 3: simple p-distance
#'   (default).
#' @param slope Numeric, relative gap width. Only used when
#'   `method = "abgd"`. Default to 1.5.
#' @param webserver A result file obtained from the ABGD (`.txt`) or ASAP
#'   (`.csv`) web server. When given, no executable is required. Default to
#'   NULL.
#' @param outfolder Path to the folder where the external program writes its
#'   output. Default to NULL, i.e. a temporary location.
#' @param align Logical, if TRUE (default) reference sequences of unequal
#'   length are aligned with [align_pq()] before being submitted.
#'   Both ABGD and ASAP require an aligned FASTA and silently return nothing
#'   otherwise, so set `align = FALSE` only when `refseq` is already aligned.
#' @param align_method Aligner used when `align = TRUE`, either `"decipher"`
#'   (default, pure R) or `"mafft"` (external program, much faster on large
#'   `refseq` slots). See [align_pq()].
#' @param mafft_exec Path to the MAFFT executable. Only used when
#'   `align_method = "mafft"`. Default to NULL, i.e. the usual lookup of
#'   [is_mafft_installed()].
#' @param merge_taxa Logical, if TRUE (default) taxa of the same partition are
#'   merged and a phyloseq object is returned. If FALSE, the partition table is
#'   returned untouched.
#' @param tax_adjust Handling of taxonomic disagreements within a partition.
#'   See [MiscMetabar::merge_taxa_vec()]. Default to 1.
#' @param rank_propagation Logical, default TRUE, whether bad ranks are
#'   propagated to lower ranks. See [MiscMetabar::merge_taxa_vec()].
#' @param keep_temporary_files Logical, if TRUE the temporary FASTA file
#'   submitted to the external program is kept. Default to FALSE.
#' @param verbose Logical, if TRUE report the number of taxa before and after
#'   delimitation. Default to FALSE.
#'
#' @return A \code{\link[phyloseq]{phyloseq-class}} object whose taxa are the
#'   delimited species, or a data.frame with columns `taxa` and `partition`
#'   when `merge_taxa = FALSE`.
#'
#' @details ABGD and ASAP are not available on Windows. Taxa that the external
#'   program does not return a partition for are dropped by
#'   [MiscMetabar::merge_taxa_vec()] and reported with a warning.
#'
#'   Both programs need an aligned FASTA. Metabarcoding reference sequences are
#'   almost never aligned, hence the `align = TRUE` default; alignment of a
#'   large `refseq` slot can be slow with the default `"decipher"` backend, in
#'   which case `align_method = "mafft"` is usually much faster.
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
#' @seealso [delim_multi_pq()], [is_delim_installed()], [align_pq()],
#'   [MiscMetabar::postcluster_pq()], [MiscMetabar::merge_taxa_vec()]
#'
#' @examplesIf phylopq::is_delim_installed("asap")
#' library(MiscMetabar)
#' data(data_fungi_mini)
#' pq_asap <- delim_pq(data_fungi_mini, method = "asap", verbose = TRUE, slope=0.5)
#' phyloseq::ntaxa(pq_asap)
#'
#' \dontrun{
#' # Give an explicit path to the executable and inspect the partitions only
#' partitions <- delim_pq(
#'   data_fungi_mini,
#'   method = "abgd",
#'   exe = "/usr/local/bin/abgd",
#'   slope = 0.5,
#'   merge_taxa = FALSE
#' )
#' head(partitions)
#'
#' # Reuse a result file downloaded from the ASAP web server
#' pq <- delim_pq(data_fungi_mini, method = "asap", webserver = "asap.csv")
#'
#' # Align with MAFFT instead of DECIPHER, much faster on large refseq slots
#' pq <- delim_pq(data_fungi_mini, method = "asap", align_method = "mafft")
#' }
#' @importFrom MiscMetabar verify_pq merge_taxa_vec
#' @export
delim_pq <- function(
  physeq,
  method = c("asap", "abgd"),
  exe = NULL,
  model = 3,
  slope = 1.5,
  webserver = NULL,
  outfolder = NULL,
  align = TRUE,
  align_method = c("decipher", "mafft"),
  mafft_exec = NULL,
  merge_taxa = TRUE,
  tax_adjust = 1L,
  rank_propagation = TRUE,
  keep_temporary_files = FALSE,
  verbose = FALSE
) {
  verify_pq(physeq)
  method <- match.arg(method)
  align_method <- match.arg(align_method)

  if (!requireNamespace("delimtools", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg delimtools} is required by {.fn delim_pq}."
    )
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg Biostrings} is required by {.fn delim_pq}."
    )
  }
  if (is.null(physeq@refseq)) {
    cli::cli_abort(c(
      "The {.arg physeq} object has no {.field refseq} slot.",
      "i" = "{.fn delim_pq} delimits species from reference sequences."
    ))
  }

  if (is.null(webserver)) {
    exe <- resolve_delim_exe(method, exe)
  }

  dna <- align_refseq(
    Biostrings::DNAStringSet(physeq@refseq),
    align,
    align_method,
    mafft_exec
  )

  fasta_file <- file.path(tempdir(), paste0("delim_pq_", method, ".fasta"))
  Biostrings::writeXStringSet(dna, fasta_file)
  if (!keep_temporary_files) {
    on.exit(unlink(fasta_file), add = TRUE)
  }

  # Both programs write their partitions next to each other and both readers
  # pick a file up by name. Sharing one directory between runs means a run that
  # outputs nothing silently inherits the previous run's partitions, so give
  # each run its own folder unless the caller chose one.
  if (is.null(outfolder)) {
    outfolder <- tempfile("phylopq_delim_")
    dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)
    if (!keep_temporary_files) {
      on.exit(unlink(outfolder, recursive = TRUE), add = TRUE)
    }
  }

  res <- if (method == "abgd") {
    delimtools::abgd_tbl(
      infile = fasta_file,
      exe = exe,
      slope = slope,
      model = model,
      outfolder = outfolder,
      webserver = webserver,
      delimname = "abgd"
    )
  } else {
    asap_tbl_safely(
      infile = fasta_file,
      exe = exe,
      model = model,
      outfolder = outfolder,
      webserver = webserver,
      delimname = "asap",
      verbose = verbose
    )
  }

  if (is.null(res)) {
    cli::cli_abort(c(
      "{.field {method}} returned no partition.",
      "i" = "See the {.pkg delimtools} messages above for the reason."
    ))
  }

  partitions <- delim_partitions(res, physeq, method)

  if (!merge_taxa) {
    return(partitions)
  }

  # A single partition means no barcode gap was found. Merging would collapse
  # every taxon into one, which `merge_taxa_vec()` cannot do; report the real
  # cause rather than letting it fail on a missing column.
  n_partitions <- length(unique(stats::na.omit(partitions$partition)))
  if (n_partitions <= 1 && phyloseq::ntaxa(physeq) > 1) {
    cli::cli_abort(c(
      "{.field {method}} placed all {phyloseq::ntaxa(physeq)} taxa in a
       single partition.",
      "i" = "No barcode gap was found, so there is nothing to merge.",
      if (method == "abgd") {
        c("i" = "Rerun with a lower {.arg slope} (currently {slope}).")
      } else {
        c(
          "i" = "Try another {.arg model}, or run the analysis on a
                 more variable marker."
        )
      },
      "i" = "Use {.code merge_taxa = FALSE} to inspect the partition table."
    ))
  }

  groups <- partitions$partition
  names(groups) <- partitions$taxa
  n_missing <- sum(is.na(groups))
  if (n_missing > 0) {
    cli::cli_warn(
      "{n_missing} taxa {?was/were} not delimited by {.field {method}} and
       {?is/are} dropped."
    )
  }

  new_obj <- merge_taxa_vec(
    physeq,
    groups,
    tax_adjust = tax_adjust,
    rank_propagation = rank_propagation
  )

  if (verbose) {
    cli::cli_inform(c(
      "v" = "{phyloseq::ntaxa(physeq)} taxa delimited into
             {phyloseq::ntaxa(new_obj)} species by {.field {method}}."
    ))
  }

  return(new_obj)
}

#' Is an ABGD or ASAP executable available?
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Check whether the external program needed by [delim_pq()] is
#' installed, either at an explicit path or through the usual lookup. Useful to
#' guard examples, tests and vignette chunks.
#'
#' @param method One of `"asap"` (default) or `"abgd"`.
#' @param path Optional path to the executable. Default to NULL, in which case
#'   `method` is looked up in three places, in order: the `phylopq.asappath` /
#'   `phylopq.abgdpath` option, a copy installed by [install_asap()] /
#'   [install_abgd()], then the system `PATH`.
#'
#' @return A logical of length one. FALSE when the `delimtools` package is not
#'   installed, so that the check also guards the R-level dependency.
#'
#' @author Adrien Taudière
#'
#' @seealso [delim_pq()]
#'
#' @examples
#' is_delim_installed("asap")
#' is_delim_installed("abgd")
#' @export
is_delim_installed <- function(method = c("asap", "abgd"), path = NULL) {
  method <- match.arg(method)
  if (!requireNamespace("delimtools", quietly = TRUE)) {
    return(FALSE)
  }
  if (!is.null(path)) {
    return(file.exists(path))
  }
  exe <- find_delim_exe(method)
  nzchar(exe) && file.exists(exe)
}

#' Align reference sequences before a delimitation run
#'
#' @inheritParams delim_pq
#' @param dna A `DNAStringSet` of reference sequences.
#' @return An aligned `DNAStringSet`.
#' @noRd
#' @keywords internal
align_refseq <- function(
  dna,
  align,
  align_method = "decipher",
  mafft_exec = NULL
) {
  if (length(unique(Biostrings::width(dna))) == 1) {
    return(dna)
  }
  if (!align) {
    cli::cli_abort(c(
      "The {.field refseq} sequences have different lengths.",
      "i" = "ABGD and ASAP require an aligned FASTA.",
      "i" = "Use {.code align = TRUE} to align them first."
    ))
  }
  align_pq(dna, method = align_method, exec = mafft_exec)
}

#' Find the executable required by [delim_pq()]
#'
#' @inheritParams delim_pq
#' @return A length-one character path to the executable.
#' @noRd
#' @keywords internal
resolve_delim_exe <- function(method, exe) {
  if (is.null(exe)) {
    exe <- find_delim_exe(method)
  }
  if (!nzchar(exe) || !file.exists(exe)) {
    cli::cli_abort(c(
      "The {.field {method}} executable was not found.",
      "i" = "Install it with {.code phylopq::install_asap()} or
             {.code phylopq::install_abgd()}.",
      "i" = "Give its path with {.arg exe}, or a web-server result file with
             {.arg webserver}.",
      "i" = "Installation instructions:
             {.url https://bioinfo.mnhn.fr/abi/public/{method}/}."
    ))
  }
  exe
}

#' Run ASAP and read the partition file it actually wrote
#'
#' `delimtools::asap_tbl()` cannot drive current ASAP builds, for two reasons:
#'
#' * it looks for `{outfolder}/{basename(infile)}.Partition_1.csv`, but ASAP
#'   derives that prefix itself and may shorten it depending on the input, so
#'   the file is missed. `asap_tbl()` then *silently* falls back to a table
#'   placing every sequence in one partition, which would merge all taxa into
#'   a single species without any warning;
#' * it finishes by tidying two "rogue" files that older ASAP builds dropped in
#'   the working directory, resolving them with [tools::file_path_as_absolute()]
#'   — which errors when they are absent — before testing `file.exists()`.
#'   Current ASAP writes everything to its `-o` folder, so a successful run
#'   aborts with `file '<name>.res.cvs' does not exist`.
#'
#' So the external program is invoked here, from a scratch working directory,
#' and the partition file is located by pattern rather than by guessed name.
#' Parsing is still delegated to `asap_tbl()` through its `webserver` argument,
#' which returns early and therefore avoids both problems. A web-server result
#' supplied by the user is passed straight through.
#'
#' @inheritParams delim_pq
#' @param infile Path to the FASTA file passed to ASAP.
#' @param delimname Name given to the delimitation column.
#' @return A `tbl_df` with columns `labels` and `delimname`.
#' @noRd
#' @keywords internal
asap_tbl_safely <- function(
  infile,
  exe,
  model,
  outfolder,
  webserver,
  delimname,
  verbose = FALSE
) {
  # A web-server result needs no local run.
  if (!is.null(webserver)) {
    return(delimtools::asap_tbl(
      infile = infile,
      exe = exe,
      model = model,
      outfolder = outfolder,
      webserver = webserver,
      delimname = delimname
    ))
  }

  # Resolve while the original working directory is still current.
  infile <- normalizePath(infile, mustWork = TRUE)
  if (is.null(outfolder)) {
    outfolder <- tempdir()
  }
  outfolder <- normalizePath(outfolder, mustWork = TRUE)

  scratch <- file.path(tempdir(), "phylopq_asap_cwd")
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  old_wd <- setwd(scratch)
  on.exit(setwd(old_wd), add = TRUE)

  asap_log <- suppressWarnings(system2(
    exe,
    args = c(
      "-d",
      model,
      "-a",
      "-o",
      shQuote(outfolder),
      shQuote(infile)
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  if (verbose) {
    cli::cli_verbatim(asap_log)
  }
  status <- attr(asap_log, "status")
  if (!is.null(status) && status != 0) {
    cli::cli_abort(c(
      "{.field asap} failed with status {status}.",
      "x" = paste(utils::tail(asap_log, 10), collapse = "\n")
    ))
  }

  # ASAP names its output from a prefix it derives itself, so match the
  # pattern and take the most recent hit rather than trusting a built name.
  hits <- list.files(
    outfolder,
    pattern = "\\.Partition_1\\.csv$",
    full.names = TRUE
  )
  if (length(hits) == 0) {
    cli::cli_abort(c(
      "{.field asap} produced no partition file in {.path {outfolder}}.",
      "x" = paste(utils::tail(asap_log, 10), collapse = "\n")
    ))
  }
  partition <- hits[order(file.mtime(hits), decreasing = TRUE)][[1]]

  res <- delimtools::asap_tbl(
    infile = infile,
    webserver = partition,
    delimname = delimname
  )
  # The webserver branch reads every column as character.
  res[[delimname]] <- as.integer(res[[delimname]])
  res
}

#' Map a delimtools partition table onto the taxa of a phyloseq object
#'
#' @inheritParams delim_pq
#' @param res A `tbl_df` returned by `delimtools::abgd_tbl()` or
#'   `delimtools::asap_tbl()`.
#' @return A data.frame with columns `taxa` and `partition`, ordered as
#'   `taxa_names(physeq)`.
#' @noRd
#' @keywords internal
delim_partitions <- function(res, physeq, method) {
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  if (!"labels" %in% colnames(res)) {
    cli::cli_abort(
      "The table returned by {.pkg delimtools} has no {.field labels} column."
    )
  }
  partition_col <- setdiff(colnames(res), "labels")
  if (length(partition_col) == 0) {
    cli::cli_abort(
      "The table returned by {.pkg delimtools} has no partition column."
    )
  }

  taxa <- phyloseq::taxa_names(physeq)
  data.frame(
    taxa = taxa,
    partition = res[[partition_col[1]]][match(taxa, res$labels)],
    stringsAsFactors = FALSE
  )
}
