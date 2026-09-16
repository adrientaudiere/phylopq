################################################################################
#' Git repositories the placement programs are built from
#'
#' Both EPA-ng and gappa are C++ projects built on the `genesis` library, which
#' they carry as a git **submodule**. A GitHub release tarball does not contain
#' submodules, so — unlike ABGD and ASAP, which are plain C tarballs — these two
#' can only be fetched with a recursive `git clone`.
#'
#' @param tool One of `"epa-ng"` or `"gappa"`.
#' @return A length-one character URL.
#' @noRd
#' @keywords internal
place_source_repo <- function(tool = c("epa-ng", "gappa")) {
  tool <- match.arg(tool)
  switch(
    tool,
    "epa-ng" = "https://github.com/pierrebarbera/epa-ng.git",
    "gappa" = "https://github.com/lczech/gappa.git"
  )
}

################################################################################
#' Check that the tools needed to build a placement program are present
#'
#' EPA-ng bundles `pll-modules`, which in turn bundles `libpll`, whose CMake
#' build calls `BISON_TARGET()` and `FLEX_TARGET()`. Missing `bison` or `flex`
#' therefore fails at the *configure* step with an unhelpful
#' `Unknown CMake command "BISON_TARGET"`, several minutes into the clone, so
#' they are checked up front. gappa builds on `genesis` alone, which carries no
#' parser generator, and needs neither.
#'
#' @inheritParams place_source_repo
#' @return NULL, invisibly, or an error listing what is missing.
#' @noRd
#' @keywords internal
check_place_build_tools <- function(tool = c("epa-ng", "gappa")) {
  tool <- match.arg(tool)

  needed <- c("git", "cmake", "make")
  if (tool == "epa-ng") {
    needed <- c(needed, "bison", "flex")
  }
  missing_tools <- needed[!nzchar(Sys.which(needed))]

  compiler <- c("g++", "c++", "clang++")
  if (!any(nzchar(Sys.which(compiler)))) {
    missing_tools <- c(missing_tools, "a C++ compiler (g++ or clang++)")
  }

  if (length(missing_tools) > 0) {
    apt <- if (tool == "epa-ng") {
      "sudo apt install build-essential cmake git bison flex"
    } else {
      "sudo apt install build-essential cmake git"
    }
    brew <- if (tool == "epa-ng") {
      "brew install cmake bison flex"
    } else {
      "brew install cmake"
    }
    cli::cli_abort(c(
      "Cannot build {.field {tool}}: {missing_tools} not found.",
      if (tool == "epa-ng" && any(c("bison", "flex") %in% missing_tools)) {
        c(
          "i" = "{.field epa-ng} bundles {.field libpll}, whose CMake build
                 needs {.field bison} and {.field flex}; without them the
                 build fails with
                 {.code Unknown CMake command \"BISON_TARGET\"}."
        )
      },
      "i" = "On Debian/Ubuntu install them with {.code {apt}}.",
      "i" = "On macOS (Homebrew): {.code {brew}}.",
      "i" = "Or skip the build entirely with
             {.code conda install -c bioconda {tool}}."
    ))
  }
  invisible(NULL)
}

################################################################################
#' Clone, build and install EPA-ng or gappa
#'
#' @inheritParams place_source_repo
#' @param path Directory in which the `bin/` sub-directory is created.
#' @param src Optional local source directory to build instead of cloning.
#' @param force If TRUE, rebuild even when the executable is already present.
#' @param verbose If TRUE, let git, cmake and make write to the console.
#' @return The path to the installed executable (invisibly).
#' @noRd
#' @keywords internal
install_place_software <- function(
  tool = c("epa-ng", "gappa"),
  path = tools::R_user_dir("phylopq", "data"),
  src = NULL,
  force = FALSE,
  verbose = FALSE
) {
  tool <- match.arg(tool)
  # `epa-ng` holds a dash, which is neither a valid option name nor part of the
  # installer's name; both drop it.
  bare <- sub("-", "", tool)
  installer <- paste0("install_", bare)
  opt_name <- paste0("phylopq.", bare, "path")

  dest_bin_dir <- delim_bin_dir(path)
  dest_bin <- file.path(dest_bin_dir, tool)

  if (file.exists(dest_bin) && !force) {
    cli::cli_inform(c(
      "v" = "{.field {tool}} is already installed at {.path {dest_bin}}.",
      "i" = "Use {.code force = TRUE} to reinstall."
    ))
    return(invisible(dest_bin))
  }

  if (.Platform$OS.type == "windows") {
    cli::cli_abort(c(
      "{.field {tool}} cannot be installed automatically on Windows.",
      "i" = "Build it under WSL, or install it with
             {.code conda install -c bioconda {tool}}.",
      "i" = "Then point phylopq at it with
             {.code options({opt_name} = \"C:/path/to/{tool}.exe\")}."
    ))
  }

  check_place_build_tools(tool)

  # 1. Assemble the sources: a copy of a local directory, or a recursive clone.
  build_root <- file.path(tempdir(), paste0(sub("-", "_", tool), "_build"))
  unlink(build_root, recursive = TRUE)
  dir.create(build_root, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(build_root, recursive = TRUE), add = TRUE)

  if (!is.null(src) && dir.exists(src)) {
    file.copy(
      list.files(src, full.names = TRUE, all.files = TRUE, no.. = TRUE),
      build_root,
      recursive = TRUE
    )
  } else {
    repo <- if (is.null(src)) place_source_repo(tool) else src
    cli::cli_inform(c(
      "i" = "Cloning {.url {repo}} (recursively, this pulls the genesis
             submodule and takes a few minutes)."
    ))
    status <- system2(
      "git",
      c(
        "clone",
        "--recursive",
        "--depth",
        "1",
        shQuote(repo),
        shQuote(build_root)
      ),
      stdout = if (verbose) "" else FALSE,
      stderr = if (verbose) "" else FALSE
    )
    if (status != 0) {
      cli::cli_abort(c(
        "Could not clone {.field {tool}} from {.url {repo}}.",
        "i" = "Clone it by hand, then pass its directory with {.arg src},
               e.g. {.code {installer}(src = \"~/src/{tool}\")}.",
        "i" = "Or install it with {.code conda install -c bioconda {tool}}."
      ))
    }
  }

  # 2. Build. Both projects ship a top-level Makefile that drives cmake.
  cli::cli_inform(c(
    "i" = "Compiling {.field {tool}}; this usually takes several minutes."
  ))
  build_log <- suppressWarnings(system2(
    "make",
    args = c("-C", shQuote(build_root)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(build_log, "status")
  if (!is.null(status) && status != 0) {
    cli::cli_abort(c(
      "Compilation of {.field {tool}} failed with status {status}.",
      "x" = paste(utils::tail(build_log, 10), collapse = "\n"),
      "i" = "Install it with {.code conda install -c bioconda {tool}}
             instead."
    ))
  }

  # 3. Both write the executable into `bin/` of the source tree.
  built <- list.files(
    build_root,
    pattern = paste0("^", tool, "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  built <- built[!file.info(built)$isdir]
  if (length(built) == 0) {
    cli::cli_abort(
      "Compilation succeeded but no {.field {tool}} executable was produced."
    )
  }

  dir.create(dest_bin_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(built[[1]], dest_bin, overwrite = TRUE)
  Sys.chmod(dest_bin, "0755")

  cli::cli_inform(c(
    "v" = "{.field {tool}} installed at {.path {dest_bin}}."
  ))
  invisible(dest_bin)
}

################################################################################
#' Install the EPA-ng phylogenetic placement program
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Clone, compile and install EPA-ng (Barbera et al. 2019) into the phylopq
#' user data directory, where [is_epang_installed()] and [place_pq()] look for
#' it. No change to your `PATH` is needed.
#'
#' EPA-ng is a C++ project that carries its libraries as git submodules, so it
#' is fetched with a recursive `git clone` rather than a release tarball, and
#' built with cmake. That needs `git`, `cmake`, `make`, a C++ compiler, and —
#' because the bundled `libpll` generates a parser — `bison` and `flex`. On
#' Debian/Ubuntu:
#' `sudo apt install build-essential cmake git bison flex`. The build takes
#' several minutes.
#'
#' **`conda install -c bioconda epa-ng` is quicker** whenever conda is
#' available; this function exists for machines where it is not.
#'
#' @param path (default: `tools::R_user_dir("phylopq", "data")`) Directory in
#'   which the `bin/` sub-directory holding the executable is created.
#' @param src (default: NULL) Optional local source directory to build instead
#'   of cloning, or an alternative git URL. A directory is copied rather than
#'   built in place.
#' @param force (default: FALSE) If TRUE, rebuild even when EPA-ng is already
#'   installed.
#' @param verbose (default: FALSE) If TRUE, let git, cmake and make write to
#'   the console.
#'
#' @return The path to the installed executable, invisibly.
#'
#' @author Adrien Taudière
#'
#' @references
#' Barbera P., Kozlov A.M., Czech L., Morel B., Darriba D., Flouri T.,
#' Stamatakis A. (2019) EPA-ng: Massively Parallel Evolutionary Placement of
#' Genetic Sequences. *Systematic Biology* 68(2):365-369.
#' \doi{10.1093/sysbio/syy054}
#'
#' @seealso [place_pq()], [is_epang_installed()], [install_gappa()]
#'
#' @examples
#' \dontrun{
#' install_epang()
#' is_epang_installed()
#'
#' # Rebuild from a clone of your own
#' install_epang(src = "~/src/epa-ng", force = TRUE, verbose = TRUE)
#' }
#' @export
install_epang <- function(
  path = tools::R_user_dir("phylopq", "data"),
  src = NULL,
  force = FALSE,
  verbose = FALSE
) {
  install_place_software(
    "epa-ng",
    path = path,
    src = src,
    force = force,
    verbose = verbose
  )
}

################################################################################
#' Install the gappa placement-analysis program
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Clone, compile and install gappa (Czech et al. 2020) into the phylopq user
#' data directory, where [is_gappa_installed()] and [assign_placement_pq()]
#' look for it. No change to your `PATH` is needed.
#'
#' gappa is a C++ project that carries the `genesis` library as a git
#' submodule, so it is fetched with a recursive `git clone` rather than a
#' release tarball, and built with cmake. That needs `git`, `cmake`, `make` and
#' a C++ compiler — but not the `bison` and `flex` that [install_epang()]
#' additionally requires. The build takes several minutes.
#'
#' **`conda install -c bioconda gappa` is quicker** whenever conda is
#' available; this function exists for machines where it is not.
#'
#' @inheritParams install_epang
#'
#' @return The path to the installed executable, invisibly.
#'
#' @author Adrien Taudière
#'
#' @references
#' Czech L., Barbera P., Stamatakis A. (2020) Genesis and Gappa: processing,
#' analyzing and visualizing phylogenetic (placement) data. *Bioinformatics*
#' 36(10):3263-3265. \doi{10.1093/bioinformatics/btaa070}
#'
#' @seealso [assign_placement_pq()], [is_gappa_installed()], [install_epang()]
#'
#' @examples
#' \dontrun{
#' install_gappa()
#' is_gappa_installed()
#'
#' # Rebuild from a clone of your own
#' install_gappa(src = "~/src/gappa", force = TRUE, verbose = TRUE)
#' }
#' @export
install_gappa <- function(
  path = tools::R_user_dir("phylopq", "data"),
  src = NULL,
  force = FALSE,
  verbose = FALSE
) {
  install_place_software(
    "gappa",
    path = path,
    src = src,
    force = force,
    verbose = verbose
  )
}
