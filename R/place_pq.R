#' Place the reference sequences of a phyloseq object on a reference tree
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Run a *phylogenetic placement* of the `refseq` slot of a
#' \code{\link[phyloseq]{phyloseq-class}} object onto a fixed reference
#' phylogeny with EPA-ng (Barbera et al. 2019). Rather than building a tree
#' from the metabarcoding reads themselves — which short, fast-evolving markers
#' rarely support — each query sequence is grafted onto the branch of a curated
#' reference tree it fits best, with a likelihood weight attached.
#'
#' The three steps are chained for you:
#'
#' 1. the query sequences are aligned **into** the reference alignment with
#'    `mafft --add --keeplength`, so that they end up on exactly the reference's
#'    columns, which EPA-ng requires;
#' 2. EPA-ng places them on the reference tree and writes a `.jplace` file;
#' 3. the `.jplace` file is read back with [BoSSA::read_jplace()] and flattened
#'    into a data.frame.
#'
#' Use [ref_package_pq()] to assemble the reference and [assign_placement_pq()]
#' to turn the placements into a taxonomy.
#'
#' @section The reference is the science:
#' Placement never questions the reference tree — it only asks where a query
#' fits on it. Everything therefore rests on that reference, and the two halves
#' of it answer to different constraints:
#'
#' * **The reference alignment must match your primers.** It has to be the same
#'   locus, over a comparable region, as the amplicon you sequenced. This is
#'   the constraint that cannot be relaxed: placing ITS2 reads on an ITS1
#'   reference yields confident-looking nonsense.
#' * **The reference tree need not come from that locus at all.** Being fixed,
#'   it is the natural place to inject evidence a single barcode cannot carry —
#'   a multi-locus phylogeny, a topology constrained by morphology, ecology or
#'   mating data, the published tree of an integrative-taxonomy study. Short
#'   markers rarely resolve deep relationships; a tree built from several loci
#'   and morpho-ecological characters does, and the placement inherits that
#'   resolution.
#'
#' In practice the reference package is almost always **external** to the
#' phyloseq object: curated vouchers, not your own ASVs. Three arrangements are
#' supported, all through [ref_package_pq()]:
#'
#' 1. a fully external reference — sequences *and* tree from a database or a
#'    published study, every taxon of `physeq` a query;
#' 2. external reference sequences with the tree inferred from them, when no
#'    published phylogeny exists for the group;
#' 3. a mixed reference — a few securely identified taxa of your own dataset
#'    added to an external reference to thicken the clade your queries fall
#'    into, the rest placed with `exclude_taxa`.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   with a non-empty `refseq` slot, or a `DNAStringSet` of query sequences.
#' @param ref_alignment The reference alignment, either a path to an aligned
#'   FASTA file or a `DNAStringSet` whose sequences all share a width. It must
#'   cover the **same marker as your queries**. Required unless `refpkg` is
#'   given.
#' @param ref_tree The reference tree, either a path to a Newick file or a
#'   `phylo` object. Its tip labels must match the names of `ref_alignment`.
#'   It is taken as fixed and never re-estimated, so it may — and ideally
#'   should — come from more evidence than the barcode alone. Required unless
#'   `refpkg` is given.
#' @param refpkg A `phylopq_refpkg` object from [ref_package_pq()], carrying
#'   the reference alignment, tree and model together. Give either this or
#'   `ref_alignment` + `ref_tree`.
#' @param query_taxa Character vector of taxa of `physeq` to place. Default to
#'   NULL, i.e. all of them.
#' @param exclude_taxa Character vector of taxa of `physeq` **not** to place.
#'   Use it for the taxa you put into the reference package itself, which would
#'   otherwise be placed on a branch they define.
#' @param model The evolutionary model handed to EPA-ng, either a model string
#'   such as `"GTR+G"` (default) or the path to a RAxML-ng `.bestModel` file.
#'   A model fitted on the reference tree is strongly preferable to the
#'   default.
#' @param query Optional path to a FASTA file of query sequences **already
#'   aligned** on the reference's columns. Default to NULL, i.e. the `refseq`
#'   slot of `physeq` is aligned for you with MAFFT.
#' @param exec Path to the EPA-ng executable. Default to NULL, in which case it
#'   is looked up with [is_epang_installed()]: first the `phylopq.epangpath`
#'   option, then a copy installed by [install_epang()], then the system
#'   `PATH`.
#' @param mafft_exec Path to the MAFFT executable, used to align the queries.
#'   Default to NULL, i.e. the `MiscMetabar.mafftpath` option then the system
#'   `PATH`. Ignored when `query` is given.
#' @param outdir Directory EPA-ng writes into. Default to NULL, i.e. a
#'   temporary location. The `.jplace` file is named `epa_result.jplace` there.
#' @param threads Number of threads given to EPA-ng and to MAFFT. Default to 1.
#' @param epang_args A length-one character vector of further arguments passed
#'   on to EPA-ng, e.g. `"--filter-acc-lwr 0.99"`. Default to `""`.
#' @param best_only Logical, if TRUE (default) only the placement with the
#'   highest likelihood weight ratio is kept for each query. If FALSE every
#'   candidate branch is returned.
#' @param return_jplace Logical, if TRUE return the `jplace` object of
#'   [BoSSA::read_jplace()] itself instead of the flattened data.frame. Default
#'   to FALSE.
#' @param keep_temporary_files Logical, if TRUE the FASTA files handed to MAFFT
#'   and EPA-ng are kept. Default to FALSE.
#' @param cmd_is_run Logical, if FALSE the commands are returned as a named
#'   character vector instead of being run. Default to TRUE. Useful to inspect,
#'   log or adapt the invocation, and to document the workflow without the
#'   external programs installed.
#' @param verbose Logical, if TRUE report each step and let EPA-ng and MAFFT
#'   write to the console. Default to FALSE.
#'
#' @return A data.frame with one row per placement and the columns `query` (the
#'   taxa name), `edge` (the branch index in the returned tree), `jplace_edge`
#'   (the branch number as EPA-ng numbered it), `like_weight_ratio`,
#'   `likelihood`, `distal_length` and `pendant_length`, carrying the path to
#'   the `.jplace` file as its `jplace_file` attribute and the reference tree as
#'   its `tree` attribute.
#'
#'   A `jplace` object when `return_jplace = TRUE`, or a named character vector
#'   of commands when `cmd_is_run = FALSE`.
#'
#' @details EPA-ng is an external program and is not installed with phylopq.
#'   The quickest route is conda (`conda install -c bioconda epa-ng`);
#'   [install_epang()] builds it from source instead.
#'
#'   The likelihood weight ratio (LWR) of a placement is the share of the total
#'   likelihood that branch accounts for; it is the number to filter on. A
#'   query spread thinly over many branches is placed with little confidence,
#'   whatever its best branch is, which is why `best_only = FALSE` is worth
#'   inspecting before trusting an assignment.
#'
#' @author Adrien Taudière
#'
#' @references
#' Barbera P., Kozlov A.M., Czech L., Morel B., Darriba D., Flouri T.,
#' Stamatakis A. (2019) EPA-ng: Massively Parallel Evolutionary Placement of
#' Genetic Sequences. *Systematic Biology* 68(2):365-369.
#' \doi{10.1093/sysbio/syy054}
#'
#' Matsen F.A., Kodner R.B., Armbrust E.V. (2010) pplacer: linear time
#' maximum-likelihood and Bayesian phylogenetic placement of sequences onto a
#' fixed reference tree. *BMC Bioinformatics* 11:538.
#' \doi{10.1186/1471-2105-11-538}
#'
#' @seealso [ref_package_pq()], [assign_placement_pq()],
#'   [is_epang_installed()], [install_epang()], [MiscMetabar::align_pq()],
#'   [BoSSA::read_jplace()]
#'
#' @examples
#' # The commands are built without EPA-ng or MAFFT being installed
#' \donttest{
#' library(MiscMetabar)
#' data(data_fungi_mini)
#' df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 12000)
#'
#' place_pq(
#'   df,
#'   ref_alignment = "reference_alignment.fasta",
#'   ref_tree = "reference_tree.newick",
#'   model = "GTR+G",
#'   cmd_is_run = FALSE
#' )
#' }
#'
#' \dontrun{
#' # Workflow 1 - a fully external reference package: curated sequences for
#' # your marker, and a tree from a multi-locus integrative study.
#' refpkg <- ref_package_pq(
#'   "unite_refs_ITS2.fasta",
#'   tree = "integrative_multilocus.newick",
#'   model = "reference.bestModel"
#' )
#' placements <- place_pq(df, refpkg = refpkg, threads = 4, verbose = TRUE)
#' head(placements)
#' attr(placements, "jplace_file")
#'
#' # Workflow 2 - external sequences, no published tree for the group
#' refpkg <- ref_package_pq("refs.fasta", tree_method = "ml")
#' placements <- place_pq(df, refpkg = refpkg)
#'
#' # Workflow 3 - strengthen the reference with a few taxa you are sure of,
#' # then place only the remaining ones
#' sure_taxa <- c("ASV1", "ASV12")
#' refpkg <- ref_package_pq(
#'   "refs.fasta",
#'   physeq = df,
#'   ref_taxa = sure_taxa
#' )
#' placements <- place_pq(df, refpkg = refpkg, exclude_taxa = sure_taxa)
#'
#' # The alignment and the tree can still be given separately
#' placements <- place_pq(
#'   df,
#'   ref_alignment = "reference_alignment.fasta",
#'   ref_tree = "reference_tree.newick"
#' )
#'
#' # Every candidate branch, to judge how confident each placement is
#' all_places <- place_pq(
#'   df,
#'   ref_alignment = "reference_alignment.fasta",
#'   ref_tree = "reference_tree.newick",
#'   best_only = FALSE
#' )
#'
#' # Discard the queries EPA-ng is not confident about
#' placements[placements$like_weight_ratio >= 0.9, ]
#'
#' # The jplace object itself, for BoSSA's plotting functions
#' jp <- place_pq(
#'   df,
#'   ref_alignment = "reference_alignment.fasta",
#'   ref_tree = "reference_tree.newick",
#'   return_jplace = TRUE
#' )
#' plot(jp, type = "precise")
#' }
#' @importFrom MiscMetabar verify_pq
#' @export
place_pq <- function(
  physeq,
  ref_alignment = NULL,
  ref_tree = NULL,
  refpkg = NULL,
  model = "GTR+G",
  query_taxa = NULL,
  exclude_taxa = NULL,
  query = NULL,
  exec = NULL,
  mafft_exec = NULL,
  outdir = NULL,
  threads = 1,
  epang_args = "",
  best_only = TRUE,
  return_jplace = FALSE,
  keep_temporary_files = FALSE,
  cmd_is_run = TRUE,
  verbose = FALSE
) {
  if (inherits(physeq, "phyloseq")) {
    verify_pq(physeq)
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg Biostrings} is required by {.fn place_pq}.")
  }

  if (!is.null(refpkg)) {
    if (!inherits(refpkg, "phylopq_refpkg")) {
      cli::cli_abort(
        "{.arg refpkg} must be a {.cls phylopq_refpkg} object as returned by
         {.fn ref_package_pq}, not {.cls {class(refpkg)}}."
      )
    }
    if (!is.null(ref_alignment) || !is.null(ref_tree)) {
      cli::cli_abort(c(
        "Give either {.arg refpkg} or {.arg ref_alignment} +
         {.arg ref_tree}, not both.",
        "i" = "{.arg refpkg} already carries the alignment and the tree."
      ))
    }
    ref_alignment <- refpkg$alignment
    ref_tree <- refpkg$tree
    if (missing(model)) {
      model <- refpkg$model
    }
  }
  if (is.null(ref_alignment) || is.null(ref_tree)) {
    cli::cli_abort(c(
      "A reference alignment and a reference tree are required.",
      "i" = "Give them with {.arg ref_alignment} and {.arg ref_tree}, or
             assemble them once with {.fn ref_package_pq} and pass the
             result as {.arg refpkg}."
    ))
  }

  physeq <- select_query_taxa(physeq, query_taxa, exclude_taxa, verbose)

  # The `.jplace` file is the result, and `assign_placement_pq()` reads it back
  # from the path carried by the returned object, so `outdir` is never deleted
  # here. A temporary one still goes away when the R session ends.
  if (is.null(outdir)) {
    outdir <- tempfile("phylopq_place_")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Only the files phylopq wrote itself are cleaned up: a path given by the
  # caller points at their own reference package and must survive the call.
  # A dry run only assembles the command, possibly for another machine, so the
  # reference files are not required to exist here.
  ref_msa_file <- as_fasta_path(
    ref_alignment,
    "phylopq_ref_msa_",
    check_exists = cmd_is_run
  )
  ref_tree_file <- as_newick_path(
    ref_tree,
    "phylopq_ref_tree_",
    check_exists = cmd_is_run
  )
  if (!keep_temporary_files) {
    written <- c(ref_msa_file, ref_tree_file)
    written <- written[vapply(
      list(ref_msa_file, ref_tree_file),
      function(x) isTRUE(attr(x, "temporary")),
      logical(1)
    )]
    if (length(written) > 0) {
      on.exit(unlink(written), add = TRUE)
    }
  }

  query_file <- query
  query_cmd <- NULL
  if (is.null(query_file)) {
    aligned <- align_query_to_ref(
      physeq,
      ref_msa_file = ref_msa_file,
      mafft_exec = mafft_exec,
      threads = threads,
      cmd_is_run = cmd_is_run,
      verbose = verbose
    )
    query_file <- aligned$path
    query_cmd <- aligned$cmd
    if (!keep_temporary_files) {
      on.exit(unlink(query_file), add = TRUE)
    }
  }

  epang_exec <- if (cmd_is_run) resolve_place_exe("epa-ng", exec) else "epa-ng"
  epang_argv <- paste(
    "--ref-msa",
    shQuote(ref_msa_file),
    "--tree",
    shQuote(ref_tree_file),
    "--query",
    shQuote(query_file),
    "--model",
    shQuote(model),
    "--outdir",
    shQuote(outdir),
    "--threads",
    threads,
    "--redo",
    epang_args
  )

  if (!cmd_is_run) {
    return(stats::setNames(
      c(query_cmd, paste(shQuote(epang_exec), epang_argv)),
      c(if (!is.null(query_cmd)) "mafft", "epa-ng")
    ))
  }

  if (verbose) {
    cli::cli_inform(c("i" = "Placing with {.field epa-ng}."))
  }
  status <- system2(
    epang_exec,
    epang_argv,
    stdout = if (verbose) "" else FALSE,
    stderr = if (verbose) "" else FALSE
  )

  jplace_file <- file.path(outdir, "epa_result.jplace")
  if (status != 0 || !file.exists(jplace_file)) {
    cli::cli_abort(c(
      "{.field epa-ng} exited with status {status} and wrote no
       {.field .jplace} file.",
      "i" = "Check that the tip labels of {.arg ref_tree} match the names of
             {.arg ref_alignment}.",
      "i" = "Use {.code verbose = TRUE} to see its output, or
             {.code cmd_is_run = FALSE} to inspect the command."
    ))
  }

  read_placements(
    jplace_file,
    best_only = best_only,
    return_jplace = return_jplace,
    verbose = verbose
  )
}

#' Assign a taxonomy from a phylogenetic placement
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Turn the `.jplace` file produced by [place_pq()] into a taxonomy, with
#' `gappa examine assign` (Czech et al. 2020). Each query is assigned the
#' consensus taxonomic path of the reference tips around the branch it was
#' placed on, so a query that falls deep inside a well-sampled clade is
#' resolved to a fine rank while one placed near the root is resolved only to a
#' coarse one.
#'
#' This is the complement of a similarity-based assignment: it answers "where
#' does this sequence sit in the reference phylogeny", not "which reference
#' sequence does it look most like".
#'
#' @param jplace (required) Path to a `.jplace` file, or a data.frame returned
#'   by [place_pq()] (its `jplace_file` attribute is used).
#' @param taxon_file (required) Path to the tab-separated file mapping each
#'   reference tip label to its taxonomic path, in the format
#'   `gappa examine assign` expects: one line per tip, the label, a tab, then
#'   the path with its ranks separated by `;`.
#' @param physeq Optional \code{\link[phyloseq]{phyloseq-class}} object. When
#'   given, the assigned taxonomy is written into its `tax_table` and the
#'   object is returned instead of the table.
#' @param ranks Character vector naming the ranks of `taxon_file`, used as the
#'   columns of the returned table. Default to the seven usual ranks.
#' @param min_alwr Numeric in \[0, 1\], the accumulated likelihood weight a
#'   clade must reach for a query to be assigned that deep. Default to 0.5.
#'   gappa reports one row per query **and per taxonomic depth**, each with the
#'   weight accumulated over that clade; the assignment kept here is the
#'   deepest row still clearing this floor, so raising `min_alwr` yields
#'   coarser but safer assignments and lowering it yields finer but shakier
#'   ones.
#' @param exec Path to the gappa executable. Default to NULL, in which case it
#'   is looked up with [is_gappa_installed()]: first the `phylopq.gappapath`
#'   option, then a copy installed by [install_gappa()], then the system
#'   `PATH`.
#' @param outdir Directory gappa writes into. Default to NULL, i.e. a temporary
#'   location.
#' @param gappa_args A length-one character vector of further arguments passed
#'   on to gappa, e.g. `"--consensus-thresh 0.9"`. Default to `""`.
#' @param cmd_is_run Logical, if FALSE the command is returned instead of being
#'   run. Default to TRUE.
#' @param verbose Logical, if TRUE report each step and let gappa write to the
#'   console. Default to FALSE.
#'
#' @return A data.frame with one row per query and the columns `query`, `alwr`
#'   (the accumulated likelihood weight of the rank it was assigned to) and one
#'   per rank of `ranks`; or the `physeq` object with that taxonomy in its
#'   `tax_table` when `physeq` is given; or the command as a length-one
#'   character vector when `cmd_is_run = FALSE`.
#'
#' @details gappa is an external program and is not installed with phylopq. The
#'   quickest route is conda (`conda install -c bioconda gappa`);
#'   [install_gappa()] builds it from source instead.
#'
#'   Queries that gappa returns no assignment for keep `NA` at every rank
#'   rather than being dropped, so the returned table always covers every query
#'   of the `.jplace` file.
#'
#'   The `--per-query-results` flag is always passed, because gappa writes
#'   `per_query.tsv` only when asked to and this function reads nothing else;
#'   `profile.tsv`, which it writes by default, aggregates over the whole
#'   sample and carries no query names.
#'
#' @author Adrien Taudière
#'
#' @references
#' Czech L., Barbera P., Stamatakis A. (2020) Genesis and Gappa: processing,
#' analyzing and visualizing phylogenetic (placement) data. *Bioinformatics*
#' 36(10):3263-3265. \doi{10.1093/bioinformatics/btaa070}
#'
#' @seealso [place_pq()], [is_gappa_installed()], [install_gappa()]
#'
#' @examples
#' \donttest{
#' # The command is built without gappa being installed
#' assign_placement_pq(
#'   jplace = "epa_result.jplace",
#'   taxon_file = "reference_taxonomy.tsv",
#'   cmd_is_run = FALSE
#' )
#' }
#'
#' \dontrun{
#' taxo <- assign_placement_pq(
#'   jplace = "epa_result.jplace",
#'   taxon_file = "reference_taxonomy.tsv",
#'   verbose = TRUE
#' )
#' head(taxo)
#'
#' # Only keep assignments the placement is confident about
#' taxo[!is.na(taxo$alwr) & taxo$alwr >= 0.9, ]
#'
#' # Coarser but safer: a clade must carry 90% of the weight to be used
#' taxo_safe <- assign_placement_pq(
#'   jplace = "epa_result.jplace",
#'   taxon_file = "reference_taxonomy.tsv",
#'   min_alwr = 0.9
#' )
#'
#' # Write the assignment straight into a phyloseq object
#' pq <- assign_placement_pq(
#'   jplace = placements,
#'   taxon_file = "reference_taxonomy.tsv",
#'   physeq = df
#' )
#' phyloseq::tax_table(pq)[1:5, ]
#'
#' # A stricter consensus
#' taxo <- assign_placement_pq(
#'   jplace = "epa_result.jplace",
#'   taxon_file = "reference_taxonomy.tsv",
#'   gappa_args = "--consensus-thresh 0.9"
#' )
#' }
#' @export
assign_placement_pq <- function(
  jplace,
  taxon_file,
  physeq = NULL,
  ranks = c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  ),
  min_alwr = 0.5,
  exec = NULL,
  outdir = NULL,
  gappa_args = "",
  cmd_is_run = TRUE,
  verbose = FALSE
) {
  jplace_file <- as_jplace_path(jplace)

  if (is.null(outdir)) {
    outdir <- tempfile("phylopq_assign_")
    on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  gappa_exec <- if (cmd_is_run) resolve_place_exe("gappa", exec) else "gappa"
  gappa_argv <- paste(
    "examine assign",
    "--jplace-path",
    shQuote(jplace_file),
    "--taxon-file",
    shQuote(taxon_file),
    "--out-dir",
    shQuote(outdir),
    "--per-query-results",
    "--allow-file-overwriting",
    gappa_args
  )

  if (!cmd_is_run) {
    return(paste(shQuote(gappa_exec), gappa_argv))
  }

  if (!file.exists(taxon_file)) {
    cli::cli_abort(c(
      "The taxon file {.path {taxon_file}} does not exist.",
      "i" = "{.arg taxon_file} maps every reference tip label to its
             taxonomic path."
    ))
  }

  if (verbose) {
    cli::cli_inform(c("i" = "Assigning with {.field gappa}."))
  }
  status <- system2(
    gappa_exec,
    gappa_argv,
    stdout = if (verbose) "" else FALSE,
    stderr = if (verbose) "" else FALSE
  )

  per_query <- file.path(outdir, "per_query.tsv")
  if (status != 0 || !file.exists(per_query)) {
    cli::cli_abort(c(
      "{.field gappa} exited with status {status} and wrote no
       {.file per_query.tsv}.",
      "i" = "Check that the labels of {.arg taxon_file} match the tip labels
             of the reference tree.",
      "i" = "Use {.code verbose = TRUE} to see its output, or
             {.code cmd_is_run = FALSE} to inspect the command."
    ))
  }

  taxo <- read_gappa_assignment(
    per_query,
    queries = jplace_queries(jplace_file),
    ranks = ranks,
    min_alwr = min_alwr
  )

  if (is.null(physeq)) {
    return(taxo)
  }

  taxo_to_pq(taxo, physeq, ranks = ranks, verbose = verbose)
}

#' Is EPA-ng available?
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Check whether the EPA-ng placement program used by [place_pq()] is
#' installed. Useful to guard examples, tests and vignette chunks.
#'
#' @param path Optional path to the EPA-ng executable. Default to NULL, in
#'   which case it is looked up in three places, in order: the
#'   `phylopq.epangpath` option, a copy installed by [install_epang()], then
#'   the system `PATH`.
#'
#' @return A logical of length one. FALSE when the `BoSSA` package is not
#'   installed, so that the check also guards the R-level dependency needed to
#'   read the result.
#'
#' @author Adrien Taudière
#'
#' @seealso [place_pq()], [install_epang()], [is_gappa_installed()]
#'
#' @examples
#' is_epang_installed()
#' @export
is_epang_installed <- function(path = NULL) {
  if (!requireNamespace("BoSSA", quietly = TRUE)) {
    return(FALSE)
  }
  if (!is.null(path)) {
    return(file.exists(path))
  }
  exe <- find_place_exe("epa-ng")
  nzchar(exe) && file.exists(exe)
}

#' Is gappa available?
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Check whether the gappa placement-analysis program used by
#' [assign_placement_pq()] is installed. Useful to guard examples, tests and
#' vignette chunks.
#'
#' @param path Optional path to the gappa executable. Default to NULL, in which
#'   case it is looked up in three places, in order: the `phylopq.gappapath`
#'   option, a copy installed by [install_gappa()], then the system `PATH`.
#'
#' @return A logical of length one.
#'
#' @author Adrien Taudière
#'
#' @seealso [assign_placement_pq()], [install_gappa()], [is_epang_installed()]
#'
#' @examples
#' is_gappa_installed()
#' @export
is_gappa_installed <- function(path = NULL) {
  if (!is.null(path)) {
    return(file.exists(path))
  }
  exe <- find_place_exe("gappa")
  nzchar(exe) && file.exists(exe)
}

#' Locate the EPA-ng or gappa executable
#'
#' Resolution order, mirroring [find_delim_exe()]:
#' 1. the `phylopq.epangpath` / `phylopq.gappapath` option,
#' 2. a copy installed by [install_epang()] / [install_gappa()],
#' 3. the system `PATH`.
#'
#' The option name drops the dash of `epa-ng`, which is not a valid name in an
#' R option list read with `getOption()`.
#'
#' @param tool One of `"epa-ng"` or `"gappa"`.
#' @return A length-one character path, or `""` when nothing was found.
#' @noRd
#' @keywords internal
find_place_exe <- function(tool = c("epa-ng", "gappa")) {
  tool <- match.arg(tool)
  opt_name <- if (tool == "epa-ng") "phylopq.epangpath" else "phylopq.gappapath"

  opt <- getOption(opt_name)
  if (!is.null(opt) && nzchar(opt)) {
    return(opt)
  }

  local_bin <- file.path(delim_bin_dir(), tool)
  if (file.exists(local_bin)) {
    return(local_bin)
  }

  unname(Sys.which(tool))
}

#' Resolve the path to the EPA-ng or gappa executable
#'
#' @inheritParams find_place_exe
#' @param exec Optional explicit path.
#' @return A length-one character path to the executable.
#' @noRd
#' @keywords internal
resolve_place_exe <- function(tool = c("epa-ng", "gappa"), exec = NULL) {
  tool <- match.arg(tool)
  if (is.null(exec)) {
    exec <- find_place_exe(tool)
  }
  if (!nzchar(exec) || !file.exists(exec)) {
    installer <- if (tool == "epa-ng") "install_epang" else "install_gappa"
    opt_name <- if (tool == "epa-ng") {
      "phylopq.epangpath"
    } else {
      "phylopq.gappapath"
    }
    cli::cli_abort(c(
      "The {.field {tool}} executable was not found.",
      "i" = "Install it with {.code conda install -c bioconda {tool}}, or
             build it from source with {.code {installer}()}.",
      "i" = "Give its path with {.arg exec}, or set it once with
             {.code options({opt_name} = \"/path/to/{tool}\")}."
    ))
  }
  exec
}

#' Restrict a query object to the taxa that should be placed
#'
#' Workflow 3 of [ref_package_pq()] puts a few securely identified taxa of the
#' dataset into the reference; those same taxa must then be kept out of the
#' queries, or they would be placed on a branch they themselves define.
#'
#' @inheritParams place_pq
#' @return `physeq`, restricted to the selected taxa.
#' @noRd
#' @keywords internal
select_query_taxa <- function(
  physeq,
  query_taxa = NULL,
  exclude_taxa = NULL,
  verbose = FALSE
) {
  if (is.null(query_taxa) && is.null(exclude_taxa)) {
    return(physeq)
  }

  all_taxa <- if (inherits(physeq, "phyloseq")) {
    phyloseq::taxa_names(physeq)
  } else {
    names(physeq)
  }

  keep <- if (is.null(query_taxa)) all_taxa else query_taxa
  unknown <- setdiff(keep, all_taxa)
  if (length(unknown) > 0) {
    cli::cli_abort(
      "{length(unknown)} of the {.arg query_taxa} {?is/are} not {?a taxon/
       taxa} of {.arg physeq}: {.val {utils::head(unknown, 5)}}."
    )
  }
  keep <- setdiff(keep, exclude_taxa)

  if (length(keep) == 0) {
    cli::cli_abort(c(
      "No taxon left to place.",
      "i" = "{.arg query_taxa} and {.arg exclude_taxa} exclude each other."
    ))
  }
  if (verbose) {
    cli::cli_inform(c(
      "i" = "{length(keep)} of the {length(all_taxa)} taxa selected as
             queries."
    ))
  }

  if (inherits(physeq, "phyloseq")) {
    phyloseq::prune_taxa(keep, physeq)
  } else {
    physeq[keep]
  }
}

#' Coerce an alignment argument to a path to a FASTA file
#'
#' @param x A path or a `DNAStringSet`.
#' @param prefix Prefix of the temporary file written for a `DNAStringSet`.
#' @param check_exists Whether a path given by the caller must already exist.
#'   FALSE for a dry run, which only assembles a command.
#' @return A length-one character path, carrying a `temporary` attribute that
#'   says whether the file was written here and may therefore be deleted.
#' @noRd
#' @keywords internal
as_fasta_path <- function(x, prefix, check_exists = TRUE) {
  if (is.character(x) && length(x) == 1) {
    if (check_exists && !file.exists(x)) {
      cli::cli_abort("The file {.path {x}} does not exist.")
    }
    attr(x, "temporary") <- FALSE
    return(x)
  }
  if (methods::is(x, "XStringSet")) {
    path <- tempfile(prefix, fileext = ".fasta")
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), path)
    attr(path, "temporary") <- TRUE
    return(path)
  }
  cli::cli_abort(
    "Expected a path to a FASTA file or a {.cls DNAStringSet},
     not {.cls {class(x)}}."
  )
}

#' Coerce a tree argument to a path to a Newick file
#'
#' @param x A path or a `phylo` object.
#' @param prefix Prefix of the temporary file written for a `phylo`.
#' @param check_exists Whether a path given by the caller must already exist.
#'   FALSE for a dry run, which only assembles a command.
#' @return A length-one character path, carrying a `temporary` attribute that
#'   says whether the file was written here and may therefore be deleted.
#' @noRd
#' @keywords internal
as_newick_path <- function(x, prefix, check_exists = TRUE) {
  if (is.character(x) && length(x) == 1) {
    if (check_exists && !file.exists(x)) {
      cli::cli_abort("The file {.path {x}} does not exist.")
    }
    attr(x, "temporary") <- FALSE
    return(x)
  }
  if (inherits(x, "phylo")) {
    path <- tempfile(prefix, fileext = ".newick")
    ape::write.tree(x, path)
    attr(path, "temporary") <- TRUE
    return(path)
  }
  cli::cli_abort(
    "Expected a path to a Newick file or a {.cls phylo} object,
     not {.cls {class(x)}}."
  )
}

#' Coerce a `jplace` argument to a path
#'
#' @param x A path, or a data.frame returned by [place_pq()].
#' @return A length-one character path.
#' @noRd
#' @keywords internal
as_jplace_path <- function(x) {
  if (inherits(x, "data.frame") && !is.null(attr(x, "jplace_file"))) {
    x <- attr(x, "jplace_file")
  }
  if (!is.character(x) || length(x) != 1) {
    cli::cli_abort(
      "{.arg jplace} must be a path to a {.field .jplace} file or the
       data.frame returned by {.fn place_pq}, not {.cls {class(x)}}."
    )
  }
  x
}

#' Align query sequences into the reference alignment
#'
#' EPA-ng requires the queries to sit on exactly the reference alignment's
#' columns, which is what `mafft --add --keeplength` produces. MAFFT returns
#' the reference and the queries together, so the reference sequences are
#' dropped again here before the file is handed to EPA-ng.
#'
#' @inheritParams place_pq
#' @param ref_msa_file Path to the reference alignment.
#' @return A list with `path`, the aligned queries, and `cmd`, the MAFFT
#'   command that produced them.
#' @noRd
#' @keywords internal
align_query_to_ref <- function(
  physeq,
  ref_msa_file,
  mafft_exec = NULL,
  threads = 1,
  cmd_is_run = TRUE,
  verbose = FALSE
) {
  query_raw <- tempfile("phylopq_query_raw_", fileext = ".fasta")
  combined <- tempfile("phylopq_query_add_", fileext = ".fasta")
  query_ali <- tempfile("phylopq_query_ali_", fileext = ".fasta")

  dna <- if (inherits(physeq, "phyloseq")) {
    if (is.null(physeq@refseq)) {
      cli::cli_abort(c(
        "The {.arg physeq} object has no {.field refseq} slot.",
        "i" = "{.fn place_pq} places reference sequences.",
        "i" = "Give an already-aligned FASTA with {.arg query} instead."
      ))
    }
    Biostrings::DNAStringSet(physeq@refseq)
  } else if (methods::is(physeq, "XStringSet")) {
    Biostrings::DNAStringSet(physeq)
  } else {
    cli::cli_abort(
      "{.arg physeq} must be a {.cls phyloseq} or {.cls DNAStringSet} object,
       not {.cls {class(physeq)}}."
    )
  }

  exec <- if (cmd_is_run) resolve_query_mafft(mafft_exec) else "mafft"
  cmd <- paste(
    shQuote(exec),
    "--add",
    shQuote(query_raw),
    "--keeplength",
    "--thread",
    threads,
    if (verbose) "" else "--quiet",
    shQuote(ref_msa_file),
    ">",
    shQuote(combined)
  )

  if (!cmd_is_run) {
    return(list(path = query_ali, cmd = cmd))
  }

  Biostrings::writeXStringSet(dna, query_raw)
  on.exit(unlink(c(query_raw, combined)), add = TRUE)

  if (verbose) {
    cli::cli_inform(c(
      "i" = "Aligning {length(dna)} quer{?y/ies} into the reference alignment
             with {.field mafft --add --keeplength}."
    ))
  }
  status <- system2(
    exec,
    paste(
      "--add",
      shQuote(query_raw),
      "--keeplength",
      "--thread",
      threads,
      if (verbose) "" else "--quiet",
      shQuote(ref_msa_file)
    ),
    stdout = combined,
    stderr = if (verbose) "" else FALSE
  )

  if (status != 0 || !file.exists(combined) || file.size(combined) == 0) {
    cli::cli_abort(c(
      "{.field mafft --add} exited with status {status} and produced no
       alignment.",
      "i" = "Check that {.path {ref_msa_file}} is an aligned FASTA.",
      "i" = "Use {.code verbose = TRUE} to see its output."
    ))
  }

  all_seqs <- Biostrings::readDNAStringSet(combined)
  queries <- all_seqs[names(all_seqs) %in% names(dna)]
  if (length(queries) != length(dna)) {
    cli::cli_warn(
      "{length(dna) - length(queries)} of the {length(dna)} quer{?y/ies}
       {?was/were} lost by {.field mafft --add}."
    )
  }
  Biostrings::writeXStringSet(queries, query_ali)

  list(path = query_ali, cmd = cmd)
}

#' Resolve the MAFFT executable used to align the queries
#'
#' [MiscMetabar::align_pq()] drives MAFFT through `ips::mafft()`, which offers
#' no `--add`, so the program is invoked directly here and its path resolved
#' the same way MiscMetabar resolves it.
#'
#' @inheritParams place_pq
#' @return A length-one character path to the executable.
#' @noRd
#' @keywords internal
resolve_query_mafft <- function(mafft_exec = NULL) {
  if (is.null(mafft_exec)) {
    opt <- getOption("MiscMetabar.mafftpath")
    mafft_exec <- if (!is.null(opt) && nzchar(opt)) {
      opt
    } else {
      unname(Sys.which("mafft"))
    }
  }
  if (!nzchar(mafft_exec) || !file.exists(mafft_exec)) {
    cli::cli_abort(c(
      "The {.field mafft} executable was not found.",
      "i" = "Install it, e.g. {.code sudo apt install mafft} or
             {.code conda install -c bioconda mafft}.",
      "i" = "Give its path with {.arg mafft_exec}, or set it once with
             {.code options(MiscMetabar.mafftpath = \"/path/to/mafft\")}.",
      "i" = "Give an already-aligned FASTA with {.arg query} to skip this
             step entirely."
    ))
  }
  mafft_exec
}

#' Read a `.jplace` file into a flat data.frame
#'
#' [BoSSA::read_jplace()] renumbers the branches to the edge indices of the
#' `phylo` it builds, and keeps the correspondence with EPA-ng's own numbering
#' in `edge_key`. Both are returned here: `edge` indexes the returned tree,
#' `jplace_edge` is what EPA-ng wrote.
#'
#' @inheritParams place_pq
#' @param jplace_file Path to a `.jplace` file.
#' @return A data.frame of placements, or a `jplace` object.
#' @noRd
#' @keywords internal
read_placements <- function(
  jplace_file,
  best_only = TRUE,
  return_jplace = FALSE,
  verbose = FALSE
) {
  if (!requireNamespace("BoSSA", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg BoSSA} is required to read a {.field .jplace} file.",
      "i" = "Install it with {.code install.packages(\"BoSSA\")}."
    ))
  }

  jp <- BoSSA::read_jplace(jplace_file)
  if (return_jplace) {
    return(jp)
  }

  pos <- jp$placement_positions
  names_by_id <- stats::setNames(
    as.character(jp$multiclass$name),
    as.character(jp$multiclass$placement_id)
  )

  # `edge_key` is a 2-row matrix: the phylo edge index over EPA-ng's own
  # branch number.
  jplace_edge <- jp$edge_key[2, ][match(pos$edge_num, jp$edge_key[1, ])]

  res <- data.frame(
    query = unname(names_by_id[as.character(pos$placement_id)]),
    edge = pos$edge_num,
    jplace_edge = jplace_edge,
    like_weight_ratio = pos$like_weight_ratio,
    likelihood = pos$likelihood,
    distal_length = pos$distal_length,
    pendant_length = pos$pendant_length,
    stringsAsFactors = FALSE
  )

  if (best_only) {
    res <- res[order(res$query, -res$like_weight_ratio), ]
    res <- res[!duplicated(res$query), ]
  }
  res <- res[order(res$query), ]
  rownames(res) <- NULL

  if (verbose) {
    cli::cli_inform(c(
      "v" = "{length(unique(res$query))} quer{?y/ies} placed on
             {length(unique(res$edge))} branch{?es}."
    ))
  }

  attr(res, "jplace_file") <- jplace_file
  attr(res, "tree") <- jp$arbre
  res
}

#' Names of the queries of a `.jplace` file
#'
#' @param jplace_file Path to a `.jplace` file.
#' @return A character vector of query names, or `character(0)`.
#' @noRd
#' @keywords internal
jplace_queries <- function(jplace_file) {
  if (!requireNamespace("BoSSA", quietly = TRUE) || !file.exists(jplace_file)) {
    return(character(0))
  }
  jp <- try(BoSSA::read_jplace(jplace_file), silent = TRUE)
  if (inherits(jp, "try-error")) {
    return(character(0))
  }
  as.character(jp$multiclass$name)
}

#' Read the `per_query.tsv` of `gappa examine assign`
#'
#' gappa writes **one row per query and per taxonomic depth**, not one row per
#' query: a query resolved to an order yields a row for its kingdom, one for
#' its phylum, and so on down. Each row carries `LWR` (the likelihood weight
#' placed exactly at that rank) and `aLWR` (the weight accumulated over that
#' clade and everything below it), so the assignment of a query is the
#' *deepest* row whose `aLWR` still clears `min_alwr`.
#'
#' The column names have changed across gappa releases, so they are matched
#' case-insensitively rather than by position, and an unexpected layout is
#' reported with the columns that were actually found. `aLWR` is absent from
#' older releases, in which case every row is eligible and the deepest one
#' wins.
#'
#' @param path Path to `per_query.tsv`.
#' @param queries Character vector of every query of the `.jplace` file, so
#'   that queries gappa left unassigned keep a row of `NA`.
#' @inheritParams assign_placement_pq
#' @return A data.frame with the columns `query`, `alwr` and one per rank.
#' @noRd
#' @keywords internal
read_gappa_assignment <- function(path, queries, ranks, min_alwr = 0.5) {
  tab <- utils::read.table(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  cols <- tolower(names(tab))
  name_col <- which(cols %in% c("name", "query"))[1]
  path_col <- which(cols %in% c("taxopath", "taxpath", "taxon", "taxonomy"))[1]
  alwr_col <- which(cols == "alwr")[1]

  if (is.na(name_col) || is.na(path_col)) {
    cli::cli_abort(c(
      "Could not find the query and taxonomic-path columns in
       {.path {path}}.",
      "i" = "Columns found: {.val {names(tab)}}.",
      "i" = "This layout of {.field gappa examine assign} is not supported;
             please report it."
    ))
  }

  res <- empty_assignment(queries, ranks)
  if (nrow(tab) == 0) {
    return(res)
  }

  query <- as.character(tab[[name_col]])
  taxopath <- as.character(tab[[path_col]])
  depth <- lengths(strsplit(taxopath, ";", fixed = TRUE))
  alwr <- if (is.na(alwr_col)) {
    rep(NA_real_, nrow(tab))
  } else {
    suppressWarnings(as.numeric(tab[[alwr_col]]))
  }

  # Keep the rows confident enough, then the deepest one of each query.
  eligible <- if (all(is.na(alwr))) {
    rep(TRUE, nrow(tab))
  } else {
    !is.na(alwr) & alwr >= min_alwr
  }
  if (!any(eligible)) {
    return(res)
  }

  ord <- order(query[eligible], -depth[eligible])
  idx <- which(eligible)[ord]
  best <- idx[!duplicated(query[idx])]

  taxo <- split_taxopath(taxopath[best], ranks)
  hit <- match(query[best], res$query)
  res[hit, ranks] <- as.data.frame(taxo, stringsAsFactors = FALSE)
  res$alwr[hit] <- alwr[best]
  res
}

#' A table of `NA` covering every query, in a stable order
#'
#' Queries gappa returns nothing for keep a row rather than being dropped, so
#' the result always covers the whole `.jplace` file.
#'
#' @param queries Character vector of query names.
#' @inheritParams assign_placement_pq
#' @return A data.frame with the columns `query`, `alwr` and one per rank.
#' @noRd
#' @keywords internal
empty_assignment <- function(queries, ranks) {
  queries <- sort(unique(as.character(queries)))
  res <- data.frame(
    query = queries,
    alwr = rep(NA_real_, length(queries)),
    stringsAsFactors = FALSE
  )
  for (r in ranks) {
    res[[r]] <- NA_character_
  }
  res
}

#' Split `;`-separated taxonomic paths into a fixed number of ranks
#'
#' @param taxopath Character vector of `;`-separated paths.
#' @inheritParams assign_placement_pq
#' @return A character matrix with one column per rank.
#' @noRd
#' @keywords internal
split_taxopath <- function(taxopath, ranks) {
  matrix(
    unlist(lapply(
      strsplit(taxopath, ";", fixed = TRUE),
      function(x) {
        x <- trimws(x)
        x[!nzchar(x)] <- NA_character_
        length(x) <- length(ranks)
        x
      }
    )),
    nrow = length(taxopath),
    ncol = length(ranks),
    byrow = TRUE,
    dimnames = list(NULL, ranks)
  )
}

#' Write an assignment table into the `tax_table` of a phyloseq object
#'
#' @param taxo A data.frame as returned by [read_gappa_assignment()].
#' @inheritParams assign_placement_pq
#' @return The `physeq` object with the assigned taxonomy.
#' @noRd
#' @keywords internal
taxo_to_pq <- function(taxo, physeq, ranks, verbose = FALSE) {
  verify_pq(physeq)
  taxa <- phyloseq::taxa_names(physeq)

  mat <- matrix(
    NA_character_,
    nrow = length(taxa),
    ncol = length(ranks),
    dimnames = list(taxa, ranks)
  )
  shared <- intersect(taxa, taxo$query)
  mat[shared, ] <- as.matrix(taxo[match(shared, taxo$query), ranks])

  if (verbose) {
    cli::cli_inform(c(
      "v" = "{length(shared)} of the {length(taxa)} taxa of {.arg physeq}
             {?was/were} assigned from the placement."
    ))
  }
  if (length(shared) == 0) {
    cli::cli_warn(c(
      "None of the {length(taxa)} taxa of {.arg physeq} appear in the
       placement.",
      "i" = "The queries were named {.val {utils::head(taxo$query, 3)}}."
    ))
  }

  phyloseq::tax_table(physeq) <- phyloseq::tax_table(mat)
  physeq
}
