skip_on_cran()

library(MiscMetabar)


test_that("align_pq and is_mafft_installed are re-exported from MiscMetabar", {
  expect_identical(align_pq, MiscMetabar::align_pq)
  expect_identical(is_mafft_installed, MiscMetabar::is_mafft_installed)
  expect_true("align_pq" %in% getNamespaceExports("phylopq"))
  expect_true("is_mafft_installed" %in% getNamespaceExports("phylopq"))
})

test_that("the re-exported align_pq still aligns a phyloseq object", {
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("Biostrings")
  df <- suppressMessages(subset_taxa_pq(
    data_fungi_mini,
    taxa_sums(data_fungi_mini) > 12000
  ))
  ali <- align_pq(df, method = "decipher")
  expect_s4_class(ali, "DNAStringSet")
  expect_length(unique(Biostrings::width(ali)), 1)
  expect_setequal(names(ali), phyloseq::taxa_names(df))
})
