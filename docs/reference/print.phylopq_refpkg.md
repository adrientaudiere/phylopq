# Print a reference package

Print a reference package

## Usage

``` r
# S3 method for class 'phylopq_refpkg'
print(x, ...)
```

## Arguments

- x:

  A `phylopq_refpkg` object, as returned by
  [`ref_package_pq()`](https://adrientaudiere.github.io/phylopq/reference/ref_package_pq.md).

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`ref_package_pq()`](https://adrientaudiere.github.io/phylopq/reference/ref_package_pq.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
library(MiscMetabar)
data(data_fungi_mini)
df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 9000)
#> Cleaning suppress 0 taxa (  ) and 6 sample(s) ( AD26-005-H_S10_MERGED.fastq.gz / CB8-019-H_S70_MERGED.fastq.gz / DY5-004-H_S97_MERGED.fastq.gz / N23-002-B_S130_MERGED.fastq.gz / NVABM0244-M_S137_MERGED.fastq.gz / T28-ABM602-B_S162_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 45
#> Number of filtered-out ASV 23
#> Number of kept ASV 22
#> Number of kept samples 131
ref_taxa <- phyloseq::taxa_names(df)[1:8]
refpkg <- ref_package_pq(
  phyloseq::prune_taxa(ref_taxa, df),
  tree = ape::rtree(length(ref_taxa), tip.label = ref_taxa)
)
print(refpkg)
#> <phylopq reference package>
#>   8 reference sequences aligned over 363 positions
#>   a tree of 8 tips, with branch lengths
#>   model: GTR+G
# }
```
