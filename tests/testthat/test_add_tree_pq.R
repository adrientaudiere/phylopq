library(MiscMetabar)
data(data_fungi_mini)

tree_mini <- taxo2tree(data_fungi_mini)

test_that("add_tree_pq accepts a phylo object", {
  pq <- add_tree_pq(data_fungi_mini, tree_mini)
  expect_s4_class(pq, "phyloseq")
  expect_false(is.null(pq@phy_tree))
  expect_equal(phyloseq::ntaxa(pq), phyloseq::ntaxa(data_fungi_mini))
  expect_setequal(
    phyloseq::taxa_names(pq),
    phyloseq::phy_tree(pq)$tip.label
  )
})

test_that("add_tree_pq requires a tree by default", {
  expect_error(add_tree_pq(data_fungi_mini), "is required")
})

test_that("add_tree_pq builds a taxonomic tree on request", {
  pq <- add_tree_pq(data_fungi_mini, use_taxo_to_build_tree = TRUE)
  expect_s4_class(pq, "phyloseq")
  expect_setequal(
    phyloseq::taxa_names(pq),
    phyloseq::phy_tree(pq)$tip.label
  )
  expect_setequal(
    phyloseq::phy_tree(pq)$tip.label,
    tree_mini$tip.label
  )
})

test_that("add_tree_pq refuses a tree alongside use_taxo_to_build_tree", {
  expect_error(
    add_tree_pq(data_fungi_mini, tree_mini, use_taxo_to_build_tree = TRUE),
    "must be NULL"
  )
})

test_that("add_tree_pq accepts a Newick file", {
  f <- file.path(tempdir(), "test_add_tree_pq.nwk")
  ape::write.tree(tree_mini, f)
  pq <- add_tree_pq(data_fungi_mini, f)
  expect_s4_class(pq, "phyloseq")
  unlink(f)
})

test_that("add_tree_pq prunes taxa absent from the tree on request", {
  small_tree <- ape::drop.tip(tree_mini, tree_mini$tip.label[1:10])
  pq <- add_tree_pq(data_fungi_mini, small_tree, prune_taxa = TRUE)
  expect_equal(phyloseq::ntaxa(pq), phyloseq::ntaxa(data_fungi_mini) - 10)
})

test_that("add_tree_pq errors by default on taxa absent from the tree", {
  small_tree <- ape::drop.tip(tree_mini, tree_mini$tip.label[1:10])
  expect_error(add_tree_pq(data_fungi_mini, small_tree), "absent from")
})

test_that("add_tree_pq drops tips absent from the phyloseq object", {
  big_tree <- ape::read.tree(
    text = paste0(
      "(",
      sub(";$", "", ape::write.tree(tree_mini)),
      ",not_a_taxon);"
    )
  )
  pq <- add_tree_pq(data_fungi_mini, big_tree)
  expect_false("not_a_taxon" %in% phyloseq::phy_tree(pq)$tip.label)
  expect_error(
    add_tree_pq(data_fungi_mini, big_tree, drop_tips = FALSE),
    "absent from"
  )
})

test_that("add_tree_pq computes branch lengths on request", {
  pq <- add_tree_pq(data_fungi_mini, tree_mini, compute_brlen = TRUE)
  expect_false(is.null(phyloseq::phy_tree(pq)$edge.length))
})

test_that("add_tree_pq refuses to overwrite without force", {
  pq <- add_tree_pq(data_fungi_mini, tree_mini)
  expect_error(add_tree_pq(pq, tree_mini), "already contains")
  expect_s4_class(add_tree_pq(pq, tree_mini, force = TRUE), "phyloseq")
})

test_that("add_tree_pq errors on non-overlapping labels", {
  odd_tree <- tree_mini
  odd_tree$tip.label <- paste0("zzz_", odd_tree$tip.label)
  expect_error(add_tree_pq(data_fungi_mini, odd_tree), "No tip label")
})

test_that("add_tree_pq roots on an outgroup and by midpoint", {
  pq <- add_tree_pq(
    data_fungi_mini,
    tree_mini,
    root = "outgroup",
    outgroup = tree_mini$tip.label[1]
  )
  expect_true(ape::is.rooted(phyloseq::phy_tree(pq)))
  expect_error(
    add_tree_pq(data_fungi_mini, tree_mini, root = "outgroup"),
    "must be given"
  )

  skip_if_not_installed("phangorn")
  pq_mid <- add_tree_pq(
    data_fungi_mini,
    tree_mini,
    root = "midpoint",
    compute_brlen = TRUE
  )
  expect_true(ape::is.rooted(phyloseq::phy_tree(pq_mid)))
})

test_that("add_tree_pq errors on a bad tree argument", {
  expect_error(add_tree_pq(data_fungi_mini, 42), "must be a")
  expect_error(
    add_tree_pq(data_fungi_mini, "no_such_file.nwk"),
    "does not exist"
  )
})

test_that("add_tree_pq passes ... on to taxo2tree", {
  pq <- add_tree_pq(
    data_fungi_mini,
    use_taxo_to_build_tree = TRUE,
    ranks = c("Domain", "Phylum", "Class")
  )
  expect_s4_class(pq, "phyloseq")
  expect_equal(phyloseq::ntaxa(pq), phyloseq::ntaxa(data_fungi_mini))
})

test_that("add_tree_pq keeps slots consistent for verify_pq", {
  for (r in c("none", "midpoint")) {
    pq <- add_tree_pq(
      data_fungi_mini,
      tree_mini,
      compute_brlen = TRUE,
      root = r
    )
    expect_no_error(MiscMetabar::verify_pq(pq))
    expect_equal(
      phyloseq::taxa_names(pq),
      phyloseq::phy_tree(pq)$tip.label
    )
  }
})

test_that("has_brlen sees through the phyloseq numeric(0) round-trip", {
  pq <- add_tree_pq(data_fungi_mini, tree_mini)
  expect_false(has_brlen(phyloseq::phy_tree(pq)))
  expect_true(has_brlen(phyloseq::phy_tree(
    add_tree_pq(data_fungi_mini, tree_mini, compute_brlen = TRUE)
  )))
})
