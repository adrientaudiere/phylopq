library(MiscMetabar)
data(data_fungi_mini)

pq_tree <- add_tree_pq(
  data_fungi_mini,
  use_taxo_to_build_tree = TRUE,
  compute_brlen = TRUE
)

test_that("phylo_glom_pq reduces the number of taxa", {
  glom <- phylo_glom_pq(pq_tree, h = 0.3)
  expect_s4_class(glom, "phyloseq")
  expect_lt(phyloseq::ntaxa(glom), phyloseq::ntaxa(pq_tree))
  expect_equal(
    sum(phyloseq::otu_table(glom)),
    sum(phyloseq::otu_table(pq_tree))
  )
})

test_that("phylo_glom_pq with h = 0 keeps every taxon", {
  glom <- phylo_glom_pq(pq_tree, h = 0)
  expect_equal(phyloseq::ntaxa(glom), phyloseq::ntaxa(pq_tree))
})

test_that("phylo_glom_pq returns a taxa-to-cluster map", {
  map <- phylo_glom_pq(pq_tree, h = 0.3, return_map = TRUE)
  expect_s3_class(map, "data.frame")
  expect_named(map, c("taxa", "cluster"))
  expect_setequal(map$taxa, phyloseq::taxa_names(pq_tree))
  expect_equal(
    length(unique(map$cluster)),
    phyloseq::ntaxa(phylo_glom_pq(pq_tree, h = 0.3))
  )
})

test_that("phylo_glom_pq is verbose on request", {
  expect_message(
    phylo_glom_pq(pq_tree, h = 0.3, verbose = TRUE),
    "agglomerated"
  )
})

test_that("phylo_glom_pq errors without a tree or branch lengths", {
  expect_error(phylo_glom_pq(data_fungi_mini, h = 0.3), "no .*phy_tree")
  pq_nobrlen <- add_tree_pq(data_fungi_mini, use_taxo_to_build_tree = TRUE)
  expect_error(phylo_glom_pq(pq_nobrlen, h = 0.3), "no branch length")
})

test_that("phylo_glom_pq validates h", {
  expect_error(phylo_glom_pq(pq_tree, h = -1), "non-negative")
  expect_error(phylo_glom_pq(pq_tree, h = c(0.1, 0.2)), "single")
})

test_that("phylo_glom_scan_pq returns one row per threshold", {
  scan <- phylo_glom_scan_pq(pq_tree, h_values = c(0, 0.2, 0.5, 1))
  expect_s3_class(scan, "data.frame")
  expect_equal(nrow(scan), 4)
  expect_named(scan, c("h", "n_taxa", "prop_taxa"))
  expect_equal(scan$n_taxa[1], phyloseq::ntaxa(pq_tree))
  expect_true(all(diff(scan$n_taxa) <= 0))
  expect_true(all(scan$prop_taxa <= 1))
})

test_that("phylo_glom_scan_pq agrees with phylo_glom_pq", {
  scan <- phylo_glom_scan_pq(pq_tree, h_values = 0.3)
  expect_equal(scan$n_taxa, phyloseq::ntaxa(phylo_glom_pq(pq_tree, h = 0.3)))
})
