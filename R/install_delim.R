################################################################################
#' Directory into which phylopq installs delimitation binaries
#'
#' @param path Parent directory. Default to the phylopq user data directory.
#' @return A length-one character path to the `bin` sub-directory.
#' @noRd
#' @keywords internal
delim_bin_dir <- function(path = tools::R_user_dir("phylopq", "data")) {
  file.path(path, "bin")
}

################################################################################
#' Locate the ABGD or ASAP executable
#'
#' Resolution order, mirroring [MiscMetabar::find_vsearch()]:
#' 1. the `phylopq.asappath` / `phylopq.abgdpath` option,
#' 2. a copy installed by [install_asap()] / [install_abgd()],
#' 3. the system `PATH`.
#'
#' @param method One of `"asap"` or `"abgd"`.
#' @return A length-one character path, or `""` when nothing was found.
#' @noRd
#' @keywords internal
find_delim_exe <- function(method = c("asap", "abgd")) {
  method <- match.arg(method)

  opt <- getOption(paste0("phylopq.", method, "path"))
  if (!is.null(opt) && nzchar(opt)) {
    return(opt)
  }

  local_bin <- file.path(delim_bin_dir(), method)
  if (file.exists(local_bin)) {
    return(local_bin)
  }

  unname(Sys.which(method))
}

################################################################################
#' Source archives tried when installing ABGD or ASAP
#'
#' The authors' server (`bioinfo.mnhn.fr`) has been unreachable since mid-2025.
#' For ABGD the Internet Archive holds the last official tarball, which is
#' pristine upstream C and is therefore tried first. No ASAP tarball was ever
#' archived, so ASAP can only come from the iTaxoTools mirror, whose sources
#' carry a Python shim removed by [neutralise_python_shim()]. The canonical URL
#' is kept last so these functions start working by themselves should the
#' server come back.
#'
#' @param method One of `"asap"` or `"abgd"`.
#' @return A character vector of URLs, in the order they should be tried.
#' @noRd
#' @keywords internal
delim_sources <- function(method = c("asap", "abgd")) {
  method <- match.arg(method)
  canonical <- paste0(
    "https://bioinfo.mnhn.fr/abi/public/",
    method,
    "/last.tgz"
  )

  if (identical(method, "abgd")) {
    return(c(
      paste0(
        "https://web.archive.org/web/20230526144924/",
        "https://bioinfo.mnhn.fr/abi/public/abgd/last.tgz"
      ),
      "https://github.com/iTaxoTools/ABGDpy/archive/refs/heads/main.tar.gz",
      canonical
    ))
  }

  c(
    "https://github.com/iTaxoTools/ASAPy/archive/refs/heads/main.tar.gz",
    canonical
  )
}

################################################################################
#' Neutralise the iTaxoTools Python shim before a command-line build
#'
#' The iTaxoTools mirrors ship `wrapio.h`, a header that redirects stdio to
#' Python (`#define printf _printf`, ...) and unconditionally includes
#' `Python.h`. The upstream `Makefile` still builds the command-line program,
#' and `wrapio.c` is *not* among its sources, so the header alone breaks the
#' build. Replacing it with an empty file restores plain stdio and the original
#' command-line behaviour.
#'
#' Only ever applied to a throwaway copy of the sources, never to files the
#' user supplied.
#'
#' @param src_dir Directory holding the sources about to be compiled.
#' @return TRUE when the shim was neutralised, FALSE otherwise.
#' @noRd
#' @keywords internal
neutralise_python_shim <- function(src_dir) {
  shim <- file.path(src_dir, "wrapio.h")
  if (!file.exists(shim)) {
    return(FALSE)
  }
  if (!any(grepl("Python.h", readLines(shim, warn = FALSE), fixed = TRUE))) {
    return(FALSE)
  }
  writeLines(
    "/* Neutralised by phylopq: restore plain stdio for the CLI build. */",
    shim
  )
  TRUE
}

################################################################################
#' Download and unpack one candidate source archive
#'
#' @param url URL of a `.tar.gz` / `.tgz` archive.
#' @param method One of `"asap"` or `"abgd"`, used to name the temporary file.
#' @param exdir Directory to unpack into.
#' @return TRUE on success, FALSE when the archive could not be fetched or
#'   unpacked.
#' @noRd
#' @keywords internal
fetch_delim_source <- function(url, method, exdir) {
  archive <- file.path(tempdir(), paste0(method, "_src.tgz"))
  unlink(archive)
  on.exit(unlink(archive), add = TRUE)

  ok <- tryCatch(
    {
      utils::download.file(url, archive, mode = "wb", quiet = TRUE)
      file.exists(archive) && file.size(archive) > 0
    },
    error = function(e) FALSE,
    warning = function(w) FALSE
  )
  if (!isTRUE(ok)) {
    return(FALSE)
  }

  unpacked <- tryCatch(
    {
      utils::untar(archive, exdir = exdir)
      TRUE
    },
    error = function(e) FALSE,
    warning = function(w) FALSE
  )
  # A Wayback error page unpacks to nothing, so check rather than trust.
  isTRUE(unpacked) && length(list.files(exdir)) > 0
}

################################################################################
#' Download, compile and install ABGD or ASAP
#'
#' Shared engine behind [install_asap()] and [install_abgd()]. Everything is
#' built inside a throwaway directory, so sources supplied through `src` are
#' copied rather than modified in place.
#'
#' @param method One of `"asap"` or `"abgd"`.
#' @param path Directory in which the `bin/` sub-directory is created.
#' @param src Optional URL, local archive or local source directory to use
#'   instead of the built-in sources.
#' @param force If TRUE, rebuild even when the executable is already present.
#' @return The path to the installed executable (invisibly).
#' @noRd
#' @keywords internal
install_delim_software <- function(
  method = c("asap", "abgd"),
  path = tools::R_user_dir("phylopq", "data"),
  src = NULL,
  force = FALSE
) {
  method <- match.arg(method)

  dest_bin_dir <- delim_bin_dir(path)
  dest_bin <- file.path(dest_bin_dir, method)

  if (file.exists(dest_bin) && !force) {
    cli::cli_inform(c(
      "v" = "{.field {method}} is already installed at {.path {dest_bin}}.",
      "i" = "Use {.code force = TRUE} to reinstall."
    ))
    return(invisible(dest_bin))
  }

  if (.Platform$OS.type == "windows") {
    cli::cli_abort(c(
      "{.field {method}} cannot be installed automatically on Windows.",
      "i" = "Its authors distribute C sources that must be compiled.",
      "i" = "Use the iTaxoTools standalone executable instead:
             {.url https://itaxotools.org/download.html}.",
      "i" = "Then point phylopq at it with
             {.code options(phylopq.{method}path = \"C:/path/to/{method}.exe\")}."
    ))
  }

  compiler <- c("gcc", "cc")[nzchar(Sys.which(c("gcc", "cc")))]
  missing_tools <- character()
  if (!nzchar(Sys.which("make"))) {
    missing_tools <- c(missing_tools, "make")
  }
  if (length(compiler) == 0) {
    missing_tools <- c(missing_tools, "a C compiler (gcc or cc)")
  }
  if (length(missing_tools) > 0) {
    cli::cli_abort(c(
      "Cannot compile {.field {method}}: {missing_tools} not found.",
      "i" = "On Debian/Ubuntu install them with
             {.code sudo apt install build-essential}."
    ))
  }

  # 1. Assemble the sources in a throwaway build directory: a copy of a local
  # directory, a local archive, or the first reachable URL among `src` or the
  # built-in fallback chain.
  build_root <- file.path(tempdir(), paste0(method, "_build"))
  unlink(build_root, recursive = TRUE)
  dir.create(build_root, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(build_root, recursive = TRUE), add = TRUE)

  if (!is.null(src) && dir.exists(src)) {
    file.copy(
      list.files(src, full.names = TRUE, all.files = TRUE, no.. = TRUE),
      build_root,
      recursive = TRUE
    )
  } else if (!is.null(src) && file.exists(src)) {
    utils::untar(src, exdir = build_root)
  } else {
    urls <- if (is.null(src)) delim_sources(method) else src
    got <- FALSE
    for (url in urls) {
      cli::cli_inform(c("i" = "Trying {.url {url}}."))
      got <- fetch_delim_source(url, method, build_root)
      if (isTRUE(got)) {
        break
      }
      cli::cli_inform(c("x" = "Unreachable, trying the next source."))
    }
    if (!isTRUE(got)) {
      cli::cli_abort(c(
        "Could not download {.field {method}} from any known source.",
        "i" = "Tried: {.url {urls}}.",
        "i" = "Download the sources by hand, then pass the archive or its
               directory with {.arg src}, e.g.
               {.code install_{method}(src = \"~/Downloads/{method}\")}."
      ))
    }
  }

  # 2. Locate the build directory. The mirrors nest the upstream sources under
  # src/<method>/, so prefer a Makefile sitting beside the method's own files.
  makefiles <- list.files(
    build_root,
    pattern = "^[Mm]akefile$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(makefiles) == 0) {
    cli::cli_abort(
      "No Makefile found in the {.field {method}} sources."
    )
  }
  preferred <- makefiles[grepl(method, dirname(makefiles), ignore.case = TRUE)]
  src_dir <- dirname(
    if (length(preferred) > 0) preferred[[1]] else makefiles[[1]]
  )

  if (neutralise_python_shim(src_dir)) {
    cli::cli_inform(c(
      "i" = "Removed the iTaxoTools Python shim to build the command-line
             program."
    ))
  }

  # 3. Compile. The upstream Makefiles hard-code CC=gcc, so pass whichever
  # compiler is actually present.
  cli::cli_inform(c("i" = "Compiling {.field {method}}."))
  build_log <- suppressWarnings(system2(
    "make",
    args = c("-C", shQuote(src_dir), paste0("CC=", compiler[[1]])),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(build_log, "status")
  if (!is.null(status) && status != 0) {
    cli::cli_abort(c(
      "Compilation of {.field {method}} failed with status {status}.",
      "x" = paste(utils::tail(build_log, 10), collapse = "\n")
    ))
  }

  built <- list.files(
    src_dir,
    pattern = paste0("^", method, "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(built) == 0) {
    cli::cli_abort(
      "Compilation succeeded but no {.field {method}} executable was produced."
    )
  }

  # 4. Install into the phylopq user data directory.
  dir.create(dest_bin_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(built[[1]], dest_bin, overwrite = TRUE)
  Sys.chmod(dest_bin, mode = "0755")

  cli::cli_inform(c(
    "v" = "{.field {method}} installed at {.path {dest_bin}}.",
    "i" = "{.fn delim_pq} will now find it automatically."
  ))
  invisible(dest_bin)
}

################################################################################
#' Install the ASAP species-delimitation program
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Download, compile and install ASAP (Assemble Species by Automatic
#' Partitioning, Puillandre et al. 2021) into the phylopq user data directory,
#' so that [delim_pq()] and [is_delim_installed()] find it automatically.
#'
#' ASAP is distributed as C sources, so a working toolchain (`make` and `gcc`
#' or `cc`) is required. On Debian or Ubuntu these come from `build-essential`.
#'
#' @param path (default: `tools::R_user_dir("phylopq", "data")`) Directory in
#'   which the `bin/` sub-directory holding the executable is created.
#' @param src (default: NULL) Optional URL, local `.tar.gz` archive, or local
#'   source directory to build from, bypassing the built-in sources. The
#'   sources are copied before building, never modified in place.
#' @param force (default: FALSE) If TRUE, rebuild even when ASAP is already
#'   installed.
#'
#' @return The path to the installed ASAP executable (invisibly).
#'
#' @details
#' The authors' server (`bioinfo.mnhn.fr`) has been unreachable since mid-2025
#' and no ASAP tarball was ever captured by the Internet Archive, so the
#' sources come from the [iTaxoTools mirror](https://github.com/iTaxoTools/ASAPy).
#' That copy carries a `wrapio.h` header redirecting stdio to Python, which
#' breaks the command-line build; phylopq neutralises it in a throwaway copy of
#' the sources before running `make`. The canonical URL is still tried last, so
#' this keeps working if the server returns.
#'
#' On Windows, install the standalone executable from
#' <https://itaxotools.org/download.html> and point phylopq at it with
#' `options(phylopq.asappath = "C:/path/to/asap.exe")`.
#'
#' @references Puillandre N., Brouillet S., Achaz G. (2021) ASAP: assemble
#'   species by automatic partitioning. *Molecular Ecology Resources*
#'   21:609-620. \doi{10.1111/1755-0998.13281}
#'
#' @author Adrien Taudière
#'
#' @seealso [install_abgd()], [is_delim_installed()], [delim_pq()]
#'
#' @examples
#' \dontrun{
#' install_asap()
#' is_delim_installed("asap")
#'
#' # From sources downloaded by hand
#' install_asap(src = "~/Downloads/ASAPy-main")
#' }
#' @export
install_asap <- function(
  path = tools::R_user_dir("phylopq", "data"),
  src = NULL,
  force = FALSE
) {
  install_delim_software(method = "asap", path = path, src = src, force = force)
}

################################################################################
#' Install the ABGD species-delimitation program
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Download, compile and install ABGD (Automatic Barcode Gap Discovery,
#' Puillandre et al. 2012) into the phylopq user data directory, so that
#' [delim_pq()] and [is_delim_installed()] find it automatically.
#'
#' ABGD is distributed as C sources, so a working toolchain (`make` and `gcc`
#' or `cc`) is required. On Debian or Ubuntu these come from `build-essential`.
#'
#' @param path (default: `tools::R_user_dir("phylopq", "data")`) Directory in
#'   which the `bin/` sub-directory holding the executable is created.
#' @param src (default: NULL) Optional URL, local `.tar.gz` archive, or local
#'   source directory to build from, bypassing the built-in sources. The
#'   sources are copied before building, never modified in place.
#' @param force (default: FALSE) If TRUE, rebuild even when ABGD is already
#'   installed.
#'
#' @return The path to the installed ABGD executable (invisibly).
#'
#' @details
#' The authors' server (`bioinfo.mnhn.fr`) has been unreachable since mid-2025.
#' The Internet Archive captured the last official ABGD tarball, which is
#' pristine upstream C, so it is tried first; the
#' [iTaxoTools mirror](https://github.com/iTaxoTools/ABGDpy) and the canonical
#' URL follow.
#'
#' On Windows, install the standalone executable from
#' <https://itaxotools.org/download.html> and point phylopq at it with
#' `options(phylopq.abgdpath = "C:/path/to/abgd.exe")`.
#'
#' @references Puillandre N., Lambert A., Brouillet S., Achaz G. (2012) ABGD,
#'   Automatic Barcode Gap Discovery for primary species delimitation.
#'   *Molecular Ecology* 21:1864-1877.
#'   \doi{10.1111/j.1365-294X.2011.05239.x}
#'
#' @author Adrien Taudière
#'
#' @seealso [install_asap()], [is_delim_installed()], [delim_pq()]
#'
#' @examples
#' \dontrun{
#' install_abgd()
#' is_delim_installed("abgd")
#'
#' # From sources downloaded by hand
#' install_abgd(src = "~/Downloads/Abgd")
#' }
#' @export
install_abgd <- function(
  path = tools::R_user_dir("phylopq", "data"),
  src = NULL,
  force = FALSE
) {
  install_delim_software(method = "abgd", path = path, src = src, force = force)
}
