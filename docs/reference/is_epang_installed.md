# Is EPA-ng available?

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Check whether the EPA-ng placement program used by
[`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md)
is installed. Useful to guard examples, tests and vignette chunks.

## Usage

``` r
is_epang_installed(path = NULL)
```

## Arguments

- path:

  Optional path to the EPA-ng executable. Default to NULL, in which case
  it is looked up in three places, in order: the `phylopq.epangpath`
  option, a copy installed by
  [`install_epang()`](https://adrientaudiere.github.io/phylopq/reference/install_epang.md),
  then the system `PATH`.

## Value

A logical of length one. FALSE when the `BoSSA` package is not
installed, so that the check also guards the R-level dependency needed
to read the result.

## See also

[`place_pq()`](https://adrientaudiere.github.io/phylopq/reference/place_pq.md),
[`install_epang()`](https://adrientaudiere.github.io/phylopq/reference/install_epang.md),
[`is_gappa_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_gappa_installed.md)

## Author

Adrien Taudière

## Examples

``` r
is_epang_installed()
#> [1] TRUE
```
