# Is gappa available?

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Check whether the gappa placement-analysis program used by
[`assign_placement_pq()`](https://adrientaudiere.github.io/phylopq/reference/assign_placement_pq.md)
is installed. Useful to guard examples, tests and vignette chunks.

## Usage

``` r
is_gappa_installed(path = NULL)
```

## Arguments

- path:

  Optional path to the gappa executable. Default to NULL, in which case
  it is looked up in three places, in order: the `phylopq.gappapath`
  option, a copy installed by
  [`install_gappa()`](https://adrientaudiere.github.io/phylopq/reference/install_gappa.md),
  then the system `PATH`.

## Value

A logical of length one.

## See also

[`assign_placement_pq()`](https://adrientaudiere.github.io/phylopq/reference/assign_placement_pq.md),
[`install_gappa()`](https://adrientaudiere.github.io/phylopq/reference/install_gappa.md),
[`is_epang_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_epang_installed.md)

## Author

Adrien Taudière

## Examples

``` r
is_gappa_installed()
#> [1] TRUE
```
