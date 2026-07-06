library(MiscMetabar)
data(data_fungi_mini)

ranks6 <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

test_that("taxo2tree returns a phylo with one tip per taxon by default", {
  tree <- taxo2tree(data_fungi_mini, ranks = ranks6)
  expect_s3_class(tree, "phylo")
  expect_equal(length(tree$tip.label), phyloseq::ntaxa(data_fungi_mini))
  expect_setequal(tree$tip.label, phyloseq::taxa_names(data_fungi_mini))
})

test_that("taxo2tree collapses identical paths when use_taxa_names = FALSE", {
  tree_full <- taxo2tree(data_fungi_mini, ranks = ranks6)
  tree_collapsed <- taxo2tree(
    data_fungi_mini,
    ranks = ranks6,
    use_taxa_names = FALSE
  )
  expect_s3_class(tree_collapsed, "phylo")
  expect_lte(length(tree_collapsed$tip.label), length(tree_full$tip.label))
})

test_that("taxo2tree honours internal_node_singletons", {
  tree_with <- taxo2tree(
    data_fungi_mini,
    ranks = ranks6,
    internal_node_singletons = TRUE
  )
  tree_without <- taxo2tree(
    data_fungi_mini,
    ranks = ranks6,
    internal_node_singletons = FALSE
  )
  expect_s3_class(tree_with, "phylo")
  expect_s3_class(tree_without, "phylo")
})

test_that("taxo2tree handles real NA values in a used rank", {
  pq <- data_fungi_mini
  tt <- as(phyloseq::tax_table(pq), "matrix")
  tt[1:3, "Genus"] <- NA
  phyloseq::tax_table(pq) <- phyloseq::tax_table(tt)
  tree <- taxo2tree(pq, ranks = ranks6)
  expect_s3_class(tree, "phylo")
  expect_equal(length(tree$tip.label), phyloseq::ntaxa(pq))
})

test_that("taxo2tree errors when a rank is missing from the tax_table", {
  expect_error(
    taxo2tree(data_fungi_mini, ranks = c("Domain", "NotARank")),
    "not present in the tax_table"
  )
})
