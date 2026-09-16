#' Assemble a reference package for phylogenetic placement
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Build the pair [place_pq()] needs — a reference **alignment** and a
#' reference **tree** sharing the same labels — from reference sequences that
#' normally come from outside the phyloseq object.
#'
#' A reference package is the scientific heart of a placement: queries can only
#' be resolved as finely as the reference around the branch they land on, so
#' the reference is where the curation effort belongs. See the
#' *Choosing a reference package* section below.
#'
#' @param x (required) The reference sequences: a path to a FASTA file, a
#'   `DNAStringSet`, or a \code{\link[phyloseq]{phyloseq-class}} object whose
#'   `refseq` slot is used. Already-aligned sequences are kept as they are
#'   unless `force_align = TRUE`.
#' @param physeq Optional \code{\link[phyloseq]{phyloseq-class}} object to draw
#'   **additional** reference sequences from, on top of `x`. Combined with
#'   `ref_taxa`, this covers the case where a few well-identified taxa of your
#'   own dataset strengthen an external reference before the remaining taxa are
#'   placed on it.
#' @param ref_taxa Character vector of taxa names of `physeq` to add to the
#'   reference. Default to NULL, i.e. every taxon of `physeq`. Ignored when
#'   `physeq` is NULL.
#' @param tree Optional reference tree, a `phylo` object or a path to a Newick
#'   or Nexus file. **Give it whenever you have one**: a tree inferred from
#'   several loci, or constrained by morphology and ecology, carries far more
#'   evidence than one built here from the barcode alone. When NULL (default) a
#'   tree is inferred from the alignment with `tree_method`.
#' @param tree_method Method used when `tree` is NULL: `"nj"` (default,
#'   neighbour-joining on a maximum-likelihood distance), `"upgma"`, or
#'   `"ml"` (a maximum-likelihood tree optimised from the NJ starting tree,
#'   much slower). All three go through the `phangorn` package.
#' @param model The evolutionary model recorded in the package and handed to
#'   EPA-ng by [place_pq()], either a model string such as `"GTR+G"` (default)
#'   or the path to a RAxML-ng `.bestModel` file.
#' @param align_method Aligner used to build the alignment, `"decipher"`
#'   (default, pure R) or `"mafft"`. See [MiscMetabar::align_pq()].
#' @param mafft_exec Path to the MAFFT executable. Only used when
#'   `align_method = "mafft"`.
#' @param force_align Logical, if TRUE sequences that already share a width are
#'   realigned anyway. Default to FALSE.
#' @param drop_tips Logical, if TRUE (default) tips of `tree` with no sequence
#'   in the alignment are dropped. If FALSE their presence is an error.
#' @param drop_seqs Logical, if TRUE sequences with no tip in `tree` are
#'   dropped from the alignment. Default to FALSE, in which case they are an
#'   error — a sequence missing from the tree is usually a labelling mistake,
#'   not something to silently discard.
#' @param verbose Logical, if TRUE report the size of the package and every
#'   reconciliation step. Default to FALSE.
#'
#' @return An object of class `phylopq_refpkg`, a list with the elements
#'   `alignment` (a `DNAStringSet`), `tree` (a `phylo`), `model` and `n_ref`.
#'   Pass it to [place_pq()] as `refpkg`.
#'
#' @section Choosing a reference package:
#' EPA-ng places a query by asking where on the reference tree its sequence
#' fits best, so the alignment and the tree play two different roles and answer
#' to two different constraints.
#'
#' * **The reference alignment must be the same marker as your queries** — the
#'   very locus your primers amplify, over a comparable region. Placing ITS2
#'   reads on an ITS1 reference, or a short amplicon on a reference trimmed
#'   elsewhere, produces confident-looking nonsense. This is the constraint
#'   that cannot be relaxed.
#' * **The reference tree does not have to come from that marker.** It is taken
#'   as fixed and never re-estimated, so it is the right place to inject
#'   evidence the barcode does not carry: a multi-locus phylogeny, a topology
#'   constrained by morphology, ecology or mating tests, the published tree of
#'   an integrative-taxonomy study. A single short barcode rarely resolves deep
#'   relationships; a tree built from several loci and checked against
#'   morphological characters does, and placement inherits that resolution.
#'
#' The usual arrangement is therefore an **external** reference package: tips
#' named for well-identified vouchers, the alignment restricted to your marker,
#' the topology taken from the best available phylogeny of the group. Pass that
#' tree through `tree` and only the alignment is built here.
#'
#' Tip labels and sequence names must match; use `drop_tips` and `drop_seqs` to
#' say what to do with the ones that do not.
#'
#' @section Three workflows:
#' 1. **Fully external reference.** Sequences and tree both come from outside
#'    the phyloseq object; its taxa are the queries. This is the common case.
#' 2. **External sequences, tree inferred here.** No published tree is
#'    available for the group, so one is built from the reference sequences
#'    with `tree_method`. Weaker, but honest as long as the reference sequences
#'    are well identified.
#' 3. **Mixed reference.** A few securely identified taxa of your own dataset
#'    are added to an external reference with `physeq` and `ref_taxa`, to
#'    thicken the part of the tree your queries fall into; the remaining taxa
#'    are then placed with `place_pq(exclude_taxa = ref_taxa)`.
#'
#' @author Adrien Taudière
#'
#' @seealso [place_pq()], [assign_placement_pq()],
#'   [MiscMetabar::align_pq()], [add_tree_pq()]
#'
#' @examples
#' \donttest{
#' library(MiscMetabar)
#' data(data_fungi_mini)
#' df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 9000)
#'
#' # Workflow 2: reference sequences without a published tree.
#' # Here they are taken from the same object for the sake of the example;
#' # in practice they come from a curated database.
#' ref_taxa <- phyloseq::taxa_names(df)[1:8]
#' refpkg <- ref_package_pq(
#'   phyloseq::prune_taxa(ref_taxa, df),
#'   verbose = TRUE
#' )
#' refpkg
#' }
#'
#' \dontrun{
#' # Workflow 1: a fully external reference package, the common case.
#' # The tree comes from a multi-locus integrative study, the alignment from
#' # the same marker as your primers.
#' refpkg <- ref_package_pq(
#'   "unite_refs_ITS2.fasta",
#'   tree = "integrative_multilocus.newick",
#'   model = "reference.bestModel"
#' )
#' placements <- place_pq(data_fungi_mini, refpkg = refpkg)
#'
#' # Workflow 3: strengthen an external reference with a few taxa of your own,
#' # then place the rest on it.
#' sure_taxa <- c("ASV1", "ASV12")
#' refpkg <- ref_package_pq(
#'   "unite_refs_ITS2.fasta",
#'   physeq = data_fungi_mini,
#'   ref_taxa = sure_taxa,
#'   align_method = "mafft"
#' )
#' placements <- place_pq(
#'   data_fungi_mini,
#'   refpkg = refpkg,
#'   exclude_taxa = sure_taxa
#' )
#'
#' # A slower but better tree when no published one exists
#' refpkg <- ref_package_pq("refs.fasta", tree_method = "ml")
#' }
#' @export
ref_package_pq <- function(
  x,
  physeq = NULL,
  ref_taxa = NULL,
  tree = NULL,
  tree_method = c("nj", "upgma", "ml"),
  model = "GTR+G",
  align_method = c("decipher", "mafft"),
  mafft_exec = NULL,
  force_align = FALSE,
  drop_tips = TRUE,
  drop_seqs = FALSE,
  verbose = FALSE
) {
  tree_method <- match.arg(tree_method)
  align_method <- match.arg(align_method)

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg Biostrings} is required by {.fn ref_package_pq}."
    )
  }

  dna <- ref_sequences(x, physeq, ref_taxa, verbose)

  alignment <- MiscMetabar::align_pq(
    dna,
    method = align_method,
    exec = mafft_exec,
    force = force_align,
    verbose = verbose
  )

  tree <- if (is.null(tree)) {
    infer_ref_tree(alignment, tree_method, verbose)
  } else {
    as_phylo(tree)
  }

  rec <- reconcile_ref(alignment, tree, drop_tips, drop_seqs, verbose)

  res <- list(
    alignment = rec$alignment,
    tree = rec$tree,
    model = model,
    n_ref = length(rec$alignment)
  )
  class(res) <- "phylopq_refpkg"

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Reference package of {res$n_ref} sequence{?s} over
             {unique(Biostrings::width(res$alignment))} position{?s}."
    ))
  }

  res
}

#' Print a reference package
#'
#' @param x A `phylopq_refpkg` object, as returned by [ref_package_pq()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @author Adrien Taudière
#'
#' @seealso [ref_package_pq()]
#'
#' @examples
#' \donttest{
#' library(MiscMetabar)
#' data(data_fungi_mini)
#' df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 9000)
#' ref_taxa <- phyloseq::taxa_names(df)[1:8]
#' refpkg <- ref_package_pq(
#'   phyloseq::prune_taxa(ref_taxa, df),
#'   tree = ape::rtree(length(ref_taxa), tip.label = ref_taxa)
#' )
#' print(refpkg)
#' }
#' @export
print.phylopq_refpkg <- function(x, ...) {
  # `cat()` rather than cli here: cli writes to the message stream, and a
  # print method belongs on stdout so that `capture.output()` and knitr see it.
  cat("<phylopq reference package>\n")
  cat(sprintf(
    "  %d reference sequences aligned over %s positions\n",
    x$n_ref,
    paste(unique(Biostrings::width(x$alignment)), collapse = "/")
  ))
  cat(sprintf(
    "  a tree of %d tips, %s branch lengths\n",
    ape::Ntip(x$tree),
    if (has_brlen(x$tree)) "with" else "without"
  ))
  cat(sprintf("  model: %s\n", x$model))
  invisible(x)
}

#' Gather the reference sequences of [ref_package_pq()]
#'
#' @inheritParams ref_package_pq
#' @return A `DNAStringSet`.
#' @noRd
#' @keywords internal
ref_sequences <- function(x, physeq = NULL, ref_taxa = NULL, verbose = FALSE) {
  dna <- if (is.character(x) && length(x) == 1) {
    if (!file.exists(x)) {
      cli::cli_abort("The file {.path {x}} does not exist.")
    }
    Biostrings::readDNAStringSet(x)
  } else {
    as_ref_dna(x)
  }

  if (!is.null(physeq)) {
    extra <- as_ref_dna(physeq)
    if (!is.null(ref_taxa)) {
      unknown <- setdiff(ref_taxa, names(extra))
      if (length(unknown) > 0) {
        cli::cli_abort(c(
          "{length(unknown)} of the {.arg ref_taxa} {?is/are} not {?a taxon/
           taxa} of {.arg physeq}: {.val {utils::head(unknown, 5)}}.",
          "i" = "{.arg ref_taxa} names taxa of {.arg physeq} to add to the
                 reference."
        ))
      }
      extra <- extra[ref_taxa]
    }
    if (verbose) {
      cli::cli_inform(c(
        "i" = "{length(extra)} sequence{?s} from {.arg physeq} added to the
               {length(dna)} of {.arg x}."
      ))
    }
    dna <- c(dna, extra)
  }

  duplicated_names <- unique(names(dna)[duplicated(names(dna))])
  if (length(duplicated_names) > 0) {
    cli::cli_abort(c(
      "{length(duplicated_names)} reference name{?s} {?is/are} duplicated:
       {.val {utils::head(duplicated_names, 5)}}.",
      "i" = "A reference package needs one sequence per label; the tree could
             not tell them apart."
    ))
  }

  if (length(dna) < 3) {
    cli::cli_abort(
      "A reference package needs at least three sequences, not {length(dna)}."
    )
  }

  dna
}

#' Coerce a reference-sequence argument to a `DNAStringSet`
#'
#' @param x A phyloseq object or an `XStringSet`.
#' @return A `DNAStringSet`.
#' @noRd
#' @keywords internal
as_ref_dna <- function(x) {
  if (inherits(x, "phyloseq")) {
    if (is.null(x@refseq)) {
      cli::cli_abort(c(
        "The phyloseq object has no {.field refseq} slot.",
        "i" = "{.fn ref_package_pq} builds a reference from sequences."
      ))
    }
    return(Biostrings::DNAStringSet(x@refseq))
  }
  if (methods::is(x, "XStringSet")) {
    return(Biostrings::DNAStringSet(x))
  }
  cli::cli_abort(
    "Expected a path to a FASTA file, a {.cls DNAStringSet} or a
     {.cls phyloseq} object, not {.cls {class(x)}}."
  )
}

#' Coerce a tree argument to a `phylo`
#'
#' @param x A `phylo` object or a path to a Newick or Nexus file.
#' @return A `phylo` object.
#' @noRd
#' @keywords internal
as_phylo <- function(x) {
  if (inherits(x, "phylo")) {
    return(x)
  }
  if (is.character(x) && length(x) == 1) {
    if (!file.exists(x)) {
      cli::cli_abort("The file {.path {x}} does not exist.")
    }
    tree <- try(ape::read.tree(x), silent = TRUE)
    if (inherits(tree, "try-error") || is.null(tree)) {
      tree <- try(ape::read.nexus(x), silent = TRUE)
    }
    if (inherits(tree, "try-error") || is.null(tree)) {
      cli::cli_abort(c(
        "Could not read a tree from {.path {x}}.",
        "i" = "Both {.fn ape::read.tree} and {.fn ape::read.nexus} failed."
      ))
    }
    return(tree)
  }
  cli::cli_abort(
    "{.arg tree} must be a {.cls phylo} object or a path to a Newick or Nexus
     file, not {.cls {class(x)}}."
  )
}

#' Infer a reference tree from the reference alignment
#'
#' Only used when no tree is given. A tree built from a single barcode is the
#' weakest kind of reference, which is why [ref_package_pq()] asks for one
#' whenever the caller has a better source.
#'
#' @inheritParams ref_package_pq
#' @param alignment An aligned `DNAStringSet`.
#' @return A rooted-free, fully dichotomous `phylo` with positive branch
#'   lengths.
#' @noRd
#' @keywords internal
infer_ref_tree <- function(alignment, tree_method, verbose = FALSE) {
  if (!requireNamespace("phangorn", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg phangorn} is required to infer a reference tree.",
      "i" = "Install it with {.code install.packages(\"phangorn\")}, or give
             a tree with {.arg tree}."
    ))
  }
  if (verbose) {
    cli::cli_inform(c(
      "i" = "Inferring a {.field {tree_method}} tree from the reference
             alignment; give {.arg tree} to use a better one."
    ))
  }

  phang <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- phangorn::dist.ml(phang)
  tree <- if (tree_method == "upgma") {
    phangorn::upgma(dm)
  } else {
    phangorn::NJ(dm)
  }

  if (tree_method == "ml") {
    fit <- phangorn::pml(phangorn::midpoint(tree), data = phang)
    tree <- phangorn::optim.pml(
      fit,
      model = "GTR",
      rearrangement = "NNI",
      control = phangorn::pml.control(trace = 0)
    )$tree
  }

  # EPA-ng needs a strictly dichotomous tree with usable branch lengths.
  tree <- ape::multi2di(tree)
  tree$edge.length[!is.finite(tree$edge.length) | tree$edge.length < 0] <- 1e-8
  tree
}

#' Make the reference alignment and the reference tree share their labels
#'
#' @inheritParams ref_package_pq
#' @param alignment An aligned `DNAStringSet`.
#' @return A list with the reconciled `alignment` and `tree`.
#' @noRd
#' @keywords internal
reconcile_ref <- function(
  alignment,
  tree,
  drop_tips = TRUE,
  drop_seqs = FALSE,
  verbose = FALSE
) {
  extra_tips <- setdiff(tree$tip.label, names(alignment))
  if (length(extra_tips) > 0) {
    if (!drop_tips) {
      cli::cli_abort(c(
        "{length(extra_tips)} tip{?s} of {.arg tree} {?has/have} no sequence
         in the reference alignment: {.val {utils::head(extra_tips, 5)}}.",
        "i" = "Use {.code drop_tips = TRUE} to drop them."
      ))
    }
    if (length(extra_tips) >= ape::Ntip(tree) - 2) {
      cli::cli_abort(c(
        "Only {ape::Ntip(tree) - length(extra_tips)} tip{?s} of {.arg tree}
         match a reference sequence.",
        "i" = "Tip labels and sequence names look unrelated; check that both
               come from the same reference package.",
        "i" = "Tips: {.val {utils::head(tree$tip.label, 3)}}.",
        "i" = "Sequences: {.val {utils::head(names(alignment), 3)}}."
      ))
    }
    tree <- ape::drop.tip(tree, extra_tips)
    if (verbose) {
      cli::cli_inform(c(
        "i" = "{length(extra_tips)} tip{?s} without a sequence dropped from
               the tree."
      ))
    }
  }

  extra_seqs <- setdiff(names(alignment), tree$tip.label)
  if (length(extra_seqs) > 0) {
    if (!drop_seqs) {
      cli::cli_abort(c(
        "{length(extra_seqs)} reference sequence{?s} {?is/are} missing from
         {.arg tree}: {.val {utils::head(extra_seqs, 5)}}.",
        "i" = "EPA-ng needs one tip per reference sequence.",
        "i" = "Use {.code drop_seqs = TRUE} to drop them, or add them to the
               tree."
      ))
    }
    alignment <- alignment[!names(alignment) %in% extra_seqs]
    if (verbose) {
      cli::cli_inform(c(
        "i" = "{length(extra_seqs)} sequence{?s} without a tip dropped from
               the alignment."
      ))
    }
  }

  if (ape::Ntip(tree) < 3) {
    cli::cli_abort(
      "Only {ape::Ntip(tree)} reference{?s} left after reconciliation; EPA-ng
       needs at least three."
    )
  }

  list(alignment = alignment[tree$tip.label], tree = tree)
}
