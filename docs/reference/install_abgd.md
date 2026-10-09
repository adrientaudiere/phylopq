# Install the ABGD species-delimitation program

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Download, compile and install ABGD (Automatic Barcode Gap Discovery,
Puillandre et al. 2012) into the phylopq user data directory, so that
[`delim_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_pq.md)
and
[`is_delim_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_delim_installed.md)
find it automatically.

ABGD is distributed as C sources, so a working toolchain (`make` and
`gcc` or `cc`) is required. On Debian or Ubuntu these come from
`build-essential`.

## Usage

``` r
install_abgd(
  path = tools::R_user_dir("phylopq", "data"),
  src = NULL,
  force = FALSE
)
```

## Arguments

- path:

  (default: `tools::R_user_dir("phylopq", "data")`) Directory in which
  the `bin/` sub-directory holding the executable is created.

- src:

  (default: NULL) Optional URL, local `.tar.gz` archive, or local source
  directory to build from, bypassing the built-in sources. The sources
  are copied before building, never modified in place.

- force:

  (default: FALSE) If TRUE, rebuild even when ABGD is already installed.

## Value

The path to the installed ABGD executable (invisibly).

## Details

The authors' server (`bioinfo.mnhn.fr`) has been unreachable since
mid-2025. The Internet Archive captured the last official ABGD tarball,
which is pristine upstream C, so it is tried first; the [iTaxoTools
mirror](https://github.com/iTaxoTools/ABGDpy) and the canonical URL
follow.

On Windows, install the standalone executable from
<https://itaxotools.org/download.html> and point phylopq at it with
`options(phylopq.abgdpath = "C:/path/to/abgd.exe")`.

## References

Puillandre N., Lambert A., Brouillet S., Achaz G. (2012) ABGD, Automatic
Barcode Gap Discovery for primary species delimitation. *Molecular
Ecology* 21:1864-1877.
[doi:10.1111/j.1365-294X.2011.05239.x](https://doi.org/10.1111/j.1365-294X.2011.05239.x)

## See also

[`install_asap()`](https://adrientaudiere.github.io/phylopq/reference/install_asap.md),
[`is_delim_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_delim_installed.md),
[`delim_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_pq.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
install_abgd()
is_delim_installed("abgd")

# From sources downloaded by hand
install_abgd(src = "~/Downloads/Abgd")
} # }
```
