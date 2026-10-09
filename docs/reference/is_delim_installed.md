# Is an ABGD or ASAP executable available?

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Check whether the external program needed by
[`delim_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_pq.md)
is installed, either at an explicit path or through the usual lookup.
Useful to guard examples, tests and vignette chunks.

## Usage

``` r
is_delim_installed(method = c("asap", "abgd"), path = NULL)
```

## Arguments

- method:

  One of `"asap"` (default) or `"abgd"`.

- path:

  Optional path to the executable. Default to NULL, in which case
  `method` is looked up in three places, in order: the
  `phylopq.asappath` / `phylopq.abgdpath` option, a copy installed by
  [`install_asap()`](https://adrientaudiere.github.io/phylopq/reference/install_asap.md)
  /
  [`install_abgd()`](https://adrientaudiere.github.io/phylopq/reference/install_abgd.md),
  then the system `PATH`.

## Value

A logical of length one. FALSE when the `delimtools` package is not
installed, so that the check also guards the R-level dependency.

## See also

[`delim_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_pq.md)

## Author

Adrien Taudière

## Examples

``` r
is_delim_installed("asap")
#> [1] TRUE
is_delim_installed("abgd")
#> [1] TRUE
```
