# Install the gappa placement-analysis program

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Clone, compile and install gappa (Czech et al. 2020) into the phylopq
user data directory, where
[`is_gappa_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_gappa_installed.md)
and
[`assign_placement_pq()`](https://adrientaudiere.github.io/phylopq/reference/assign_placement_pq.md)
look for it. No change to your `PATH` is needed.

gappa is a C++ project that carries the `genesis` library as a git
submodule, so it is fetched with a recursive `git clone` rather than a
release tarball, and built with cmake. That needs `git`, `cmake`, `make`
and a C++ compiler — but not the `bison` and `flex` that
[`install_epang()`](https://adrientaudiere.github.io/phylopq/reference/install_epang.md)
additionally requires. The build takes several minutes.

**`conda install -c bioconda gappa` is quicker** whenever conda is
available; this function exists for machines where it is not.

## Usage

``` r
install_gappa(
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

Czech L., Barbera P., Stamatakis A. (2020) Genesis and Gappa:
processing, analyzing and visualizing phylogenetic (placement) data.
*Bioinformatics* 36(10):3263-3265.
[doi:10.1093/bioinformatics/btaa070](https://doi.org/10.1093/bioinformatics/btaa070)

## See also

[`assign_placement_pq()`](https://adrientaudiere.github.io/phylopq/reference/assign_placement_pq.md),
[`is_gappa_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_gappa_installed.md),
[`install_epang()`](https://adrientaudiere.github.io/phylopq/reference/install_epang.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
install_gappa()
is_gappa_installed()

# Rebuild from a clone of your own
install_gappa(src = "~/src/gappa", force = TRUE, verbose = TRUE)
} # }
```
