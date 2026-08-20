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
#'   given. Default to NULL, in which case the executable is looked up on the
#'   `PATH` with [is_delim_installed()].
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
#'   length are aligned with [DECIPHER::AlignSeqs()] before being submitted.
#'   Both ABGD and ASAP require an aligned FASTA and silently return nothing
#'   otherwise, so set `align = FALSE` only when `refseq` is already aligned.
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
#'   large `refseq` slot can be slow.
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
#' @seealso [is_delim_installed()], [MiscMetabar::postcluster_pq()],
#'   [MiscMetabar::merge_taxa_vec()]
#'
#' @examplesIf phylopq::is_delim_installed("asap")
#' library(MiscMetabar)
#' data(data_fungi_mini)
#' pq_asap <- delim_pq(data_fungi_mini, method = "asap", verbose = TRUE)
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
  merge_taxa = TRUE,
  tax_adjust = 1L,
  rank_propagation = TRUE,
  keep_temporary_files = FALSE,
  verbose = FALSE
) {
  verify_pq(physeq)
  method <- match.arg(method)

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

  dna <- align_refseq(Biostrings::DNAStringSet(physeq@refseq), align)

  fasta_file <- file.path(tempdir(), paste0("delim_pq_", method, ".fasta"))
  Biostrings::writeXStringSet(dna, fasta_file)
  if (!keep_temporary_files) {
    on.exit(unlink(fasta_file), add = TRUE)
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
    delimtools::asap_tbl(
      infile = fasta_file,
      exe = exe,
      model = model,
      outfolder = outfolder,
      webserver = webserver,
      delimname = "asap"
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
#' installed, either at an explicit path or on the `PATH`. Useful to guard
#' examples, tests and vignette chunks.
#'
#' @param method One of `"asap"` (default) or `"abgd"`.
#' @param path Optional path to the executable. Default to NULL, in which case
#'   `method` is looked up on the `PATH` with [base::Sys.which()].
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
  nzchar(Sys.which(method))
}

#' Align reference sequences before a delimitation run
#'
#' @inheritParams delim_pq
#' @param dna A `DNAStringSet` of reference sequences.
#' @return An aligned `DNAStringSet`.
#' @noRd
#' @keywords internal
align_refseq <- function(dna, align) {
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
  if (!requireNamespace("DECIPHER", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg DECIPHER} is required to align the {.field refseq} slot."
    )
  }
  DECIPHER::AlignSeqs(dna, anchor = NA, verbose = FALSE)
}

#' Find the executable required by [delim_pq()]
#'
#' @inheritParams delim_pq
#' @return A length-one character path to the executable.
#' @noRd
#' @keywords internal
resolve_delim_exe <- function(method, exe) {
  if (is.null(exe)) {
    exe <- unname(Sys.which(method))
  }
  if (!nzchar(exe) || !file.exists(exe)) {
    cli::cli_abort(c(
      "The {.field {method}} executable was not found.",
      "i" = "Give its path with {.arg exe}, or a web-server result file with
             {.arg webserver}.",
      "i" = "Installation instructions:
             {.url https://bioinfo.mnhn.fr/abi/public/{method}/}."
    ))
  }
  exe
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
