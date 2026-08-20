test_that("find_delim_exe prefers the option over everything else", {
  old <- options(phylopq.asappath = "/nonexistent/fake-asap")
  on.exit(options(old))
  expect_identical(find_delim_exe("asap"), "/nonexistent/fake-asap")
})

test_that("find_delim_exe returns a single string", {
  res <- find_delim_exe("abgd")
  expect_type(res, "character")
  expect_length(res, 1)
})

test_that("find_delim_exe rejects an unknown method", {
  expect_error(find_delim_exe("gmyc"), "should be one of")
})

test_that("find_delim_exe finds a binary installed in the user data dir", {
  path <- file.path(tempdir(), "phylopq_find_test")
  on.exit(unlink(path, recursive = TRUE))
  bin <- file.path(path, "bin", "asap")
  dir.create(dirname(bin), recursive = TRUE, showWarnings = FALSE)
  file.create(bin)

  local_mocked_bindings(delim_bin_dir = function(...) file.path(path, "bin"))
  expect_identical(find_delim_exe("asap"), bin)
})

test_that("delim_sources puts the pristine ABGD tarball first", {
  abgd <- delim_sources("abgd")
  expect_match(abgd[[1]], "web\\.archive\\.org")
  expect_match(abgd[[length(abgd)]], "bioinfo\\.mnhn\\.fr")
})

test_that("delim_sources falls back to the mirror for ASAP", {
  asap <- delim_sources("asap")
  expect_match(asap[[1]], "iTaxoTools/ASAPy")
  expect_match(asap[[length(asap)]], "bioinfo\\.mnhn\\.fr")
  expect_false(any(grepl("web.archive.org", asap)))
})

test_that("delim_sources rejects an unknown method", {
  expect_error(delim_sources("gmyc"), "should be one of")
})

test_that("neutralise_python_shim empties a Python-tainted wrapio.h", {
  dir <- file.path(tempdir(), "phylopq_shim_test")
  on.exit(unlink(dir, recursive = TRUE))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  shim <- file.path(dir, "wrapio.h")
  writeLines(c("#include <Python.h>", "#define printf _printf"), shim)

  expect_true(neutralise_python_shim(dir))
  expect_false(any(grepl("Python.h", readLines(shim), fixed = TRUE)))
})

test_that("neutralise_python_shim leaves other sources alone", {
  dir <- file.path(tempdir(), "phylopq_shim_test2")
  on.exit(unlink(dir, recursive = TRUE))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  expect_false(neutralise_python_shim(dir))

  writeLines("#include <stdio.h>", file.path(dir, "wrapio.h"))
  expect_false(neutralise_python_shim(dir))
})

test_that("install_asap does not reinstall an existing executable", {
  path <- file.path(tempdir(), "phylopq_asap_test")
  on.exit(unlink(path, recursive = TRUE))
  bin <- file.path(path, "bin", "asap")
  dir.create(dirname(bin), recursive = TRUE, showWarnings = FALSE)
  file.create(bin)

  expect_message(res <- install_asap(path = path), "already installed")
  expect_identical(res, bin)
})

test_that("install_abgd does not reinstall an existing executable", {
  path <- file.path(tempdir(), "phylopq_abgd_test")
  on.exit(unlink(path, recursive = TRUE))
  bin <- file.path(path, "bin", "abgd")
  dir.create(dirname(bin), recursive = TRUE, showWarnings = FALSE)
  file.create(bin)

  expect_message(res <- install_abgd(path = path), "already installed")
  expect_identical(res, bin)
})

test_that("install_delim_software rejects an unknown method", {
  expect_error(install_delim_software("gmyc"), "should be one of")
})

test_that("install_asap builds without touching the supplied sources", {
  skip_if_not(nzchar(Sys.which("make")))
  skip_if_not(any(nzchar(Sys.which(c("gcc", "cc")))))

  src <- file.path(tempdir(), "phylopq_src_test")
  path <- file.path(tempdir(), "phylopq_build_test")
  on.exit(unlink(c(src, path), recursive = TRUE))
  dir.create(src, recursive = TRUE, showWarnings = FALSE)

  writeLines("int main(void) { return 0; }", file.path(src, "asap.c"))
  writeLines(
    c("CC= cc", "all: asap", "asap: asap.c", "\t$(CC) -o asap asap.c"),
    file.path(src, "Makefile")
  )

  bin <- suppressMessages(install_asap(path = path, src = src))
  expect_true(file.exists(bin))
  expect_identical(unname(file.access(bin, mode = 1)), 0L)
  expect_false(file.exists(file.path(src, "asap")))
})

test_that("install_asap builds sources carrying the Python shim", {
  skip_if_not(nzchar(Sys.which("make")))
  skip_if_not(any(nzchar(Sys.which(c("gcc", "cc")))))

  src <- file.path(tempdir(), "phylopq_shim_build")
  path <- file.path(tempdir(), "phylopq_shim_install")
  on.exit(unlink(c(src, path), recursive = TRUE))
  dir.create(src, recursive = TRUE, showWarnings = FALSE)

  writeLines(
    c(
      "#include <Python.h>",
      "int _printf(const char *f, ...);",
      "#define printf _printf"
    ),
    file.path(src, "wrapio.h")
  )
  writeLines(
    c(
      "#include <stdio.h>",
      "#include \"wrapio.h\"",
      "int main(void) { printf(\"ok\\n\"); return 0; }"
    ),
    file.path(src, "asap.c")
  )
  writeLines(
    c("CC= cc", "all: asap", "asap: asap.c", "\t$(CC) -o asap asap.c"),
    file.path(src, "Makefile")
  )

  bin <- suppressMessages(install_asap(path = path, src = src))
  expect_true(file.exists(bin))
  expect_identical(system2(bin, stdout = TRUE), "ok")
})

test_that("install_abgd reports every source it tried when all fail", {
  path <- file.path(tempdir(), "phylopq_fail_test")
  on.exit(unlink(path, recursive = TRUE))
  expect_error(
    suppressMessages(
      install_abgd(path = path, src = "https://invalid.example/nope.tar.gz")
    ),
    "Could not download"
  )
})

test_that("resolve_delim_exe points at the installers", {
  old <- options(phylopq.asappath = NULL)
  on.exit(options(old))
  expect_error(
    resolve_delim_exe("asap", "/nonexistent/asap"),
    "install_asap"
  )
})
