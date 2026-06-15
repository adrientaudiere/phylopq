test_that("phylopq package loads", {
  expect_true(requireNamespace("phylopq", quietly = TRUE))
  expect_true(requireNamespace("phyloseq", quietly = TRUE))
  expect_true(requireNamespace("ape", quietly = TRUE))
})
