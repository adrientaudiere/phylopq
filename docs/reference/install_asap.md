# Install the ASAP species-delimitation program

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Download, compile and install ASAP (Assemble Species by Automatic
Partitioning, Puillandre et al. 2021) into the phylopq user data
directory, so that
[`delim_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_pq.md)
and
[`is_delim_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_delim_installed.md)
find it automatically.

ASAP is distributed as C sources, so a working toolchain (`make` and
`gcc` or `cc`) is required. On Debian or Ubuntu these come from
`build-essential`.

## Usage

``` r
install_asap(
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

  (default: FALSE) If TRUE, rebuild even when ASAP is already installed.

## Value

The path to the installed ASAP executable (invisibly).

## Details

The authors' server (`bioinfo.mnhn.fr`) has been unreachable since
mid-2025 and no ASAP tarball was ever captured by the Internet Archive,
so the sources come from the [iTaxoTools
mirror](https://github.com/iTaxoTools/ASAPy). That copy carries a
`wrapio.h` header redirecting stdio to Python, which breaks the
command-line build; phylopq neutralises it in a throwaway copy of the
sources before running `make`. The canonical URL is still tried last, so
this keeps working if the server returns.

On Windows, install the standalone executable from
<https://itaxotools.org/download.html> and point phylopq at it with
`options(phylopq.asappath = "C:/path/to/asap.exe")`.

## References

Puillandre N., Brouillet S., Achaz G. (2021) ASAP: assemble species by
automatic partitioning. *Molecular Ecology Resources* 21:609-620.
[doi:10.1111/1755-0998.13281](https://doi.org/10.1111/1755-0998.13281)

## See also

[`install_abgd()`](https://adrientaudiere.github.io/phylopq/reference/install_abgd.md),
[`is_delim_installed()`](https://adrientaudiere.github.io/phylopq/reference/is_delim_installed.md),
[`delim_pq()`](https://adrientaudiere.github.io/phylopq/reference/delim_pq.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
install_asap()
is_delim_installed("asap")

# From sources downloaded by hand
install_asap(src = "~/Downloads/ASAPy-main")
} # }
```
