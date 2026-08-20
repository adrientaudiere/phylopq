skip_on_cran()

library(MiscMetabar)
data(data_fungi_mini)

df <- suppressMessages(subset_taxa_pq(
  data_fungi_mini,
  taxa_sums(data_fungi_mini) > 9000
))

test_that("delim_run_grid gives one abgd run per slope and one asap run", {
  grid <- delim_run_grid(c("abgd", "asap"), c(0.5, 1, 1.5))
  expect_equal(nrow(grid), 4)
  expect_equal(
    grid$label,
    c("abgd_slope0.5", "abgd_slope1", "abgd_slope1.5", "asap")
  )
  expect_false(any(grepl("[0-9]", grid$alias)))
  expect_equal(anyDuplicated(grid$alias), 0)
  expect_equal(nrow(delim_run_grid("asap", c(0.5, 1))), 1)
  expect_equal(nrow(delim_run_grid("abgd", c(0.5, 0.5))), 1)
  expect_error(delim_run_grid("abgd", numeric(0)), "non-empty")
  expect_error(delim_run_grid("abgd", c(1, NA)), "non-empty")
})

test_that("alias_letters is digit-free and unique beyond 26", {
  expect_equal(alias_letters(3), c("A", "B", "C"))
  expect_equal(alias_letters(28)[27:28], c("AA", "AB"))
  expect_false(any(grepl("[0-9]", alias_letters(60))))
  expect_equal(anyDuplicated(alias_letters(60)), 0)
})

test_that("as_delim_treedata carries the columns delim_autoplot reads", {
  skip_if_not_installed("treeio")
  skip_if_not_installed("tidytree")
  tree <- ape::compute.brlen(taxo2tree(df))
  td <- as_delim_treedata(tree)
  expect_s4_class(td, "treedata")
  expect_true(all(c("support", "posterior") %in% colnames(td@data)))
  expect_equal(nrow(td@data), ape::Ntip(tree) + tree$Nnode)
  expect_false(anyNA(td@data$support))
  expect_false(anyNA(td@data$posterior))
})

test_that("delim_multi_pq errors without a refseq slot", {
  no_seq <- df
  no_seq@refseq <- NULL
  expect_error(delim_multi_pq(no_seq), "refseq")
})

test_that("delim_multi_pq rejects a tree that is not a phylo", {
  skip_if_not(is_delim_installed("abgd"))
  expect_error(
    delim_multi_pq(df, methods = "abgd", slopes = 0.5, tree = "not_a_tree"),
    "phylo"
  )
})

test_that("delim_multi_pq compares several runs and joins them", {
  skip_if_not(is_delim_installed("abgd"))
  skip_if_not(is_delim_installed("asap"))
  res <- suppressWarnings(delim_multi_pq(
    df,
    slopes = c(0.5, 1.5),
    plot = FALSE
  ))
  expect_named(res, c("delim", "joined", "summary", "plot", "tree"))
  expect_length(res$delim, 3)
  expect_named(
    res$delim,
    c("abgd_slope0.5", "abgd_slope1.5", "asap")
  )
  expect_equal(
    colnames(res$joined),
    c("labels", "abgd_slope0.5", "abgd_slope1.5", "asap")
  )
  expect_false(anyNA(res$joined))
  expect_equal(nrow(res$summary), 3)
  expect_true(all(res$summary$n_partitions >= 1))
  expect_null(res$plot)
})

test_that("delim_multi_pq warns and returns no plot without a tree", {
  skip_if_not(is_delim_installed("abgd"))
  expect_warning(
    res <- delim_multi_pq(df, methods = "abgd", slopes = c(0.5, 1.5)),
    "No tree available"
  )
  expect_null(res$plot)
  expect_null(res$tree)
})

test_that("delim_multi_pq draws a figure when a tree is given", {
  skip_if_not(is_delim_installed("abgd"))
  skip_if_not_installed("ggtree")
  skip_if_not_installed("patchwork")
  res <- suppressWarnings(delim_multi_pq(
    df,
    methods = "abgd",
    slopes = c(0.5, 1, 1.5),
    tree = taxo2tree(df)
  ))
  expect_s3_class(res$plot, "patchwork")
  expect_s3_class(res$tree, "phylo")
  expect_setequal(res$tree$tip.label, res$joined$labels)
  expect_length(res$plot, 2)
  expect_equal(res$plot[[2]]$theme$axis.text.x$angle, 45)
})

test_that("delim_multi_pq leaves the run names horizontal on request", {
  skip_if_not(is_delim_installed("abgd"))
  skip_if_not_installed("ggtree")
  skip_if_not_installed("patchwork")
  res <- suppressWarnings(delim_multi_pq(
    df,
    methods = "abgd",
    slopes = c(0.5, 1.5),
    tree = taxo2tree(df),
    label_angle = 0
  ))
  expect_null(res$plot[[2]]$theme$axis.text.x$angle)
})

test_that("delim_multi_pq handles a single run", {
  skip_if_not(is_delim_installed("abgd"))
  res <- suppressWarnings(delim_multi_pq(
    df,
    methods = "abgd",
    slopes = 0.5,
    plot = FALSE
  ))
  expect_equal(colnames(res$joined), c("labels", "abgd_slope0.5"))
  expect_true(all(grepl("^sp", res$joined$abgd_slope0.5)))
})
