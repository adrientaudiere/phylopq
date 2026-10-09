# Install the EPA-ng phylogenetic placement program

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Clone, compile and install EPA-ng (Barbera et al. 2019) into the phylopq
user data directory, where
[`is_epang_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_epang_installed.md)
and
[`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md)
look for it. No change to your `PATH` is needed.

EPA-ng is a C++ project that carries its libraries as git submodules, so
it is fetched with a recursive `git clone` rather than a release
tarball, and built with cmake. That needs `git`, `cmake`, `make`, a C++
compiler, and — because the bundled `libpll` generates a parser —
`bison` and `flex`. On Debian/Ubuntu:
`sudo apt install build-essential cmake git bison flex`. The build takes
several minutes.

**`conda install -c bioconda epa-ng` is quicker** whenever conda is
available; this function exists for machines where it is not.

## Usage

``` r
install_epang(
  path = tools::R_user_dir("phylopq", "data"),
  src = NULL,
  force = FALSE,
  verbose = FALSE
)
```

## Arguments

- path:

  (default: `tools::R_user_dir("phylopq", "data")`) Directory in which
  the `bin/` sub-directory holding the executable is created.

- src:

  (default: NULL) Optional local source directory to build instead of
  cloning, or an alternative git URL. A directory is copied rather than
  built in place.

- force:

  (default: FALSE) If TRUE, rebuild even when EPA-ng is already
  installed.

- verbose:

  (default: FALSE) If TRUE, let git, cmake and make write to the
  console.

## Value

The path to the installed executable, invisibly.

## References

Barbera P., Kozlov A.M., Czech L., Morel B., Darriba D., Flouri T.,
Stamatakis A. (2019) EPA-ng: Massively Parallel Evolutionary Placement
of Genetic Sequences. *Systematic Biology* 68(2):365-369.
[doi:10.1093/sysbio/syy054](https://doi.org/10.1093/sysbio/syy054)

## See also

[`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md),
[`is_epang_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_epang_installed.md),
[`install_gappa()`](https://adrientaudiere.github.io/phylopq/reference/install_gappa.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
install_epang()
is_epang_installed()

# Rebuild from a clone of your own
install_epang(src = "~/src/epa-ng", force = TRUE, verbose = TRUE)
} # }
```
