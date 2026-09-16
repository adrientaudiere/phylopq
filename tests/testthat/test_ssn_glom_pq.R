skip_on_cran()

library(MiscMetabar)


df <- suppressMessages(subset_taxa_pq(
  data_fungi_mini,
  taxa_sums(data_fungi_mini) > 5000
))

pairs_fixture <- data.frame(
  from = c("a", "b", "d"),
  to = c("b", "c", "e"),
  id = c(0.99, 0.95, 0.98),
  stringsAsFactors = FALSE
)
taxa_fixture <- c("a", "b", "c", "d", "e", "f")

test_that("validate_ssn_id accepts a single number in [0, 1]", {
  expect_invisible(validate_ssn_id(0.97))
  expect_error(validate_ssn_id(1.2), "between 0 and 1")
  expect_error(validate_ssn_id(-0.1), "between 0 and 1")
  expect_error(validate_ssn_id(c(0.9, 0.95)), "single number")
  expect_error(validate_ssn_id(NA), "single number")
  expect_error(validate_ssn_id("0.97"), "single number")
})

test_that("ssn_clusters chains linked taxa and isolates the rest", {
  skip_if_not_installed("igraph")
  cl <- ssn_clusters(pairs_fixture, taxa = taxa_fixture, id = 0.9)
  expect_named(cl, taxa_fixture)
  # a-b-c chained through b, d-e linked, f alone
  expect_equal(unname(cl[c("a", "b", "c")]), rep(cl[["a"]], 3))
  expect_equal(cl[["d"]], cl[["e"]])
  expect_length(unique(cl), 3)
  expect_false(cl[["f"]] %in% cl[c("a", "d")])
})

test_that("ssn_clusters splits the chain at a higher threshold", {
  skip_if_not_installed("igraph")
  cl <- ssn_clusters(pairs_fixture, taxa = taxa_fixture, id = 0.97)
  # the b-c edge (0.95) is dropped, so c becomes its own cluster
  expect_equal(cl[["a"]], cl[["b"]])
  expect_false(cl[["c"]] == cl[["b"]])
  expect_length(unique(cl), 4)
})

test_that("ssn_clusters puts every taxon in its own cluster without edges", {
  skip_if_not_installed("igraph")
  no_edge <- pairs_fixture[0, , drop = FALSE]
  cl <- ssn_clusters(no_edge, taxa = taxa_fixture, id = 0.9)
  expect_length(unique(cl), length(taxa_fixture))
})

test_that("ssn_graph carries the identity and the membership", {
  skip_if_not_installed("igraph")
  cl <- ssn_clusters(pairs_fixture, taxa = taxa_fixture, id = 0.9)
  net <- ssn_graph(pairs_fixture, taxa_fixture, id = 0.9, clusters = cl)
  expect_s3_class(net, "igraph")
  expect_length(igraph::V(net), length(taxa_fixture))
  expect_length(igraph::E(net), nrow(pairs_fixture))
  expect_equal(igraph::V(net)$nsc, unname(cl[igraph::V(net)$name]))
  expect_setequal(igraph::E(net)$id, pairs_fixture$id)
})

test_that("ssn_graph rejects labels that are not taxa names", {
  skip_if_not_installed("igraph")
  expect_error(
    ssn_graph(pairs_fixture, taxa = c("a", "b"), id = 0.9),
    "not taxa names"
  )
})

test_that("read_uc_pairs keeps H records and drops N records", {
  uc <- tempfile(fileext = ".uc")
  writeLines(
    c(
      "H\t1\t24\t95.8\t+\t0\t0\t24M\ta\tb",
      "N\t*\t*\t*\t.\t*\t*\t*\td\t*",
      "H\t2\t24\t99.0\t+\t0\t0\t24M\tb\tc"
    ),
    uc
  )
  pairs <- read_uc_pairs(uc)
  expect_s3_class(pairs, "data.frame")
  expect_named(pairs, c("from", "to", "id"))
  expect_equal(nrow(pairs), 2)
  expect_equal(pairs$id, c(0.958, 0.99))
  expect_equal(pairs$from, c("a", "b"))
  unlink(uc)
})

test_that("read_uc_pairs returns an empty frame on an empty file", {
  uc <- tempfile(fileext = ".uc")
  file.create(uc)
  pairs <- read_uc_pairs(uc)
  expect_equal(nrow(pairs), 0)
  expect_named(pairs, c("from", "to", "id"))
  unlink(uc)
})

test_that("resolve_vsearch_exec reports a missing executable", {
  expect_error(
    resolve_vsearch_exec(tempfile("no_vsearch_here_")),
    "not found"
  )
})

test_that("ssn_glom_pq errors on bad input", {
  expect_error(ssn_glom_pq(df, id = 2), "between 0 and 1")

  no_refseq <- phyloseq::phyloseq(
    phyloseq::otu_table(df),
    phyloseq::tax_table(df)
  )
  expect_error(ssn_glom_pq(no_refseq), "refseq")

  one_taxon <- phyloseq::prune_taxa(phyloseq::taxa_names(df)[1], df)
  expect_error(ssn_glom_pq(one_taxon), "At least two taxa")
})

test_that("ssn_glom_pq agglomerates taxa into network sequence clusters", {
  skip_if_not(MiscMetabar::is_vsearch_installed())
  skip_if_not_installed("igraph")
  pq <- ssn_glom_pq(df, id = 0.9)
  expect_s4_class(pq, "phyloseq")
  expect_lte(phyloseq::ntaxa(pq), phyloseq::ntaxa(df))
  # Single-taxon clusters are kept, so no abundance is lost
  expect_equal(sum(phyloseq::otu_table(pq)), sum(phyloseq::otu_table(df)))
  expect_equal(phyloseq::nsamples(pq), phyloseq::nsamples(df))
})

test_that("ssn_glom_pq returns the map and the graph on demand", {
  skip_if_not(MiscMetabar::is_vsearch_installed())
  skip_if_not_installed("igraph")
  map <- ssn_glom_pq(df, id = 0.9, return_map = TRUE)
  expect_s3_class(map, "data.frame")
  expect_named(map, c("taxa", "cluster"))
  expect_setequal(map$taxa, phyloseq::taxa_names(df))

  net <- ssn_glom_pq(df, id = 0.9, return_graph = TRUE)
  expect_s3_class(net, "igraph")
  expect_length(igraph::V(net), phyloseq::ntaxa(df))
  expect_setequal(igraph::V(net)$nsc, unique(map$cluster))
})

test_that("ssn_glom_pq at id = 1 merges no more than a lower threshold", {
  skip_if_not(MiscMetabar::is_vsearch_installed())
  skip_if_not_installed("igraph")
  n_high <- phyloseq::ntaxa(ssn_glom_pq(df, id = 1))
  n_low <- phyloseq::ntaxa(ssn_glom_pq(df, id = 0.8))
  expect_gte(n_high, n_low)
})

test_that("ssn_glom_scan_pq scans thresholds monotonically", {
  skip_if_not(MiscMetabar::is_vsearch_installed())
  skip_if_not_installed("igraph")
  id_values <- seq(0.8, 1, by = 0.05)
  scan <- ssn_glom_scan_pq(df, id_values = id_values)
  expect_s3_class(scan, "data.frame")
  expect_named(scan, c("id", "n_edges", "n_taxa", "prop_taxa"))
  expect_equal(scan$id, id_values)
  expect_equal(scan$prop_taxa, scan$n_taxa / phyloseq::ntaxa(df))
  # More stringent thresholds can only split clusters further
  expect_false(is.unsorted(scan$n_taxa))
  expect_false(is.unsorted(rev(scan$n_edges)))
  expect_lte(max(scan$n_taxa), phyloseq::ntaxa(df))
})

test_that("ssn_glom_scan_pq agrees with ssn_glom_pq", {
  skip_if_not(MiscMetabar::is_vsearch_installed())
  skip_if_not_installed("igraph")
  scan <- ssn_glom_scan_pq(df, id_values = 0.9)
  expect_equal(scan$n_taxa, phyloseq::ntaxa(ssn_glom_pq(df, id = 0.9)))
})

test_that("ssn_glom_scan_pq errors on an empty id_values", {
  expect_error(ssn_glom_scan_pq(df, id_values = numeric(0)), "at least one")
  expect_error(ssn_glom_scan_pq(df, id_values = c(0.9, 2)), "between 0 and 1")
})
