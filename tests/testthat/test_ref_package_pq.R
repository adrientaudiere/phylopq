skip_on_cran()

library(MiscMetabar)


df <- suppressMessages(subset_taxa_pq(
  data_fungi_mini,
  taxa_sums(data_fungi_mini) > 9000
))
ref_taxa <- phyloseq::taxa_names(df)[1:8]
qry_taxa <- setdiff(phyloseq::taxa_names(df), ref_taxa)
ref_pq <- phyloseq::prune_taxa(ref_taxa, df)

toy_dna <- function(n = 5) {
  set.seed(1)
  seqs <- vapply(
    seq_len(n),
    function(i) paste(sample(c("A", "C", "G", "T"), 60, TRUE), collapse = ""),
    character(1)
  )
  Biostrings::DNAStringSet(stats::setNames(seqs, paste0("ref", seq_len(n))))
}

test_that("ref_sequences accepts a phyloseq, a DNAStringSet and a path", {
  skip_if_not_installed("Biostrings")
  expect_length(ref_sequences(ref_pq), length(ref_taxa))
  expect_length(ref_sequences(toy_dna(5)), 5)

  fa <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(toy_dna(5), fa)
  on.exit(unlink(fa))
  expect_length(ref_sequences(fa), 5)

  expect_error(ref_sequences(tempfile("absent_")), "does not exist")
  expect_error(ref_sequences(1:10), "DNAStringSet")
})

test_that("ref_sequences adds taxa of physeq to the reference", {
  skip_if_not_installed("Biostrings")
  # Workflow 3: external references thickened with a few taxa of one's own
  dna <- ref_sequences(toy_dna(5), physeq = df, ref_taxa = ref_taxa[1:2])
  expect_length(dna, 7)
  expect_true(all(ref_taxa[1:2] %in% names(dna)))

  # Without ref_taxa every taxon of physeq joins the reference
  expect_length(
    ref_sequences(toy_dna(5), physeq = ref_pq),
    5 + length(ref_taxa)
  )
})

test_that("ref_sequences rejects unknown taxa and duplicated labels", {
  skip_if_not_installed("Biostrings")
  expect_error(
    ref_sequences(toy_dna(5), physeq = df, ref_taxa = "not_a_taxon"),
    "not a taxon|not taxa"
  )
  dup <- c(toy_dna(5), toy_dna(5))
  expect_error(ref_sequences(dup), "duplicated")
  expect_error(ref_sequences(toy_dna(2)), "at least three")
})

test_that("as_phylo reads a phylo, a Newick file and a Nexus file", {
  tree <- ape::rtree(5)
  expect_identical(as_phylo(tree), tree)

  nwk <- tempfile(fileext = ".newick")
  ape::write.tree(tree, nwk)
  expect_equal(ape::Ntip(as_phylo(nwk)), 5)

  nex <- tempfile(fileext = ".nex")
  ape::write.nexus(tree, file = nex)
  expect_equal(ape::Ntip(as_phylo(nex)), 5)

  on.exit(unlink(c(nwk, nex)))
  expect_error(as_phylo(tempfile("absent_")), "does not exist")
  expect_error(as_phylo(1:10), "phylo")
})

test_that("reconcile_ref drops tips without a sequence", {
  ali <- toy_dna(5)
  tree <- ape::rtree(7, tip.label = c(names(ali), "ghost1", "ghost2"))

  rec <- reconcile_ref(ali, tree, drop_tips = TRUE, drop_seqs = FALSE)
  expect_equal(ape::Ntip(rec$tree), 5)
  expect_setequal(rec$tree$tip.label, names(ali))
  # The alignment comes back in tip order, as EPA-ng expects
  expect_equal(names(rec$alignment), rec$tree$tip.label)

  expect_error(
    reconcile_ref(ali, tree, drop_tips = FALSE, drop_seqs = FALSE),
    "no sequence"
  )
  expect_message(
    reconcile_ref(ali, tree, drop_tips = TRUE, verbose = TRUE),
    "dropped from the tree"
  )
})

test_that("reconcile_ref refuses to silently drop sequences", {
  ali <- toy_dna(6)
  tree <- ape::rtree(4, tip.label = names(ali)[1:4])

  expect_error(
    reconcile_ref(ali, tree, drop_seqs = FALSE),
    "missing from"
  )
  rec <- reconcile_ref(ali, tree, drop_seqs = TRUE)
  expect_length(rec$alignment, 4)
  expect_equal(ape::Ntip(rec$tree), 4)
})

test_that("reconcile_ref reports unrelated label sets", {
  ali <- toy_dna(5)
  tree <- ape::rtree(6, tip.label = paste0("other", 1:6))
  expect_error(
    reconcile_ref(ali, tree, drop_tips = TRUE),
    "look unrelated"
  )
})

test_that("infer_ref_tree builds a usable tree", {
  skip_if_not_installed("phangorn")
  skip_if_not_installed("DECIPHER")
  ali <- MiscMetabar::align_pq(ref_pq, method = "decipher", force = TRUE)

  for (m in c("nj", "upgma")) {
    tree <- infer_ref_tree(ali, m)
    expect_s3_class(tree, "phylo")
    expect_equal(ape::Ntip(tree), length(ref_taxa))
    expect_true(ape::is.binary(tree))
    expect_true(all(tree$edge.length >= 0))
  }
})

test_that("ref_package_pq builds a package from sequences alone", {
  skip_if_not_installed("phangorn")
  skip_if_not_installed("DECIPHER")
  rp <- ref_package_pq(ref_pq)

  expect_s3_class(rp, "phylopq_refpkg")
  expect_named(rp, c("alignment", "tree", "model", "n_ref"))
  expect_s4_class(rp$alignment, "DNAStringSet")
  expect_s3_class(rp$tree, "phylo")
  expect_equal(rp$n_ref, length(ref_taxa))
  expect_equal(names(rp$alignment), rp$tree$tip.label)
  expect_length(unique(Biostrings::width(rp$alignment)), 1)
  expect_equal(rp$model, "GTR+G")
})

test_that("ref_package_pq takes an external tree over inferring one", {
  skip_if_not_installed("DECIPHER")
  # Workflow 1: the tree comes from elsewhere (multi-locus, morphology, ...)
  external <- ape::rtree(length(ref_taxa), tip.label = ref_taxa)
  rp <- ref_package_pq(ref_pq, tree = external, model = "GTR+G+I")

  expect_equal(rp$model, "GTR+G+I")
  expect_setequal(rp$tree$tip.label, ref_taxa)
  # phangorn is never needed when a tree is given
  expect_equal(ape::Ntip(rp$tree), length(ref_taxa))
})

test_that("ref_package_pq mixes external sequences with taxa of physeq", {
  skip_if_not_installed("DECIPHER")
  # Workflow 3
  own <- ref_taxa[1:3]
  external <- Biostrings::DNAStringSet(phyloseq::refseq(
    phyloseq::prune_taxa(ref_taxa[4:8], df)
  ))
  names(external) <- paste0("EXT_", names(external))
  tree <- ape::rtree(8, tip.label = c(names(external), own))

  rp <- ref_package_pq(external, physeq = df, ref_taxa = own, tree = tree)
  expect_equal(rp$n_ref, 8)
  expect_true(all(own %in% names(rp$alignment)))
  expect_true(all(names(external) %in% names(rp$alignment)))
})

test_that("print.phylopq_refpkg summarises the package", {
  skip_if_not_installed("DECIPHER")
  external <- ape::rtree(length(ref_taxa), tip.label = ref_taxa)
  rp <- ref_package_pq(ref_pq, tree = external)
  expect_output(print(rp), "phylopq reference package")
  expect_output(print(rp), "GTR\\+G")
  expect_output(print(rp), "8 reference sequences")
  expect_output(print(rp), "8 tips")
  expect_invisible(print(rp))
})

test_that("select_query_taxa keeps and excludes taxa", {
  expect_equal(
    phyloseq::ntaxa(select_query_taxa(df, query_taxa = qry_taxa)),
    length(qry_taxa)
  )
  expect_equal(
    phyloseq::ntaxa(select_query_taxa(df, exclude_taxa = ref_taxa)),
    length(qry_taxa)
  )
  # No selection at all leaves the object untouched
  expect_identical(select_query_taxa(df), df)

  dna <- toy_dna(5)
  expect_length(select_query_taxa(dna, exclude_taxa = "ref1"), 4)
})

test_that("select_query_taxa errors on unknown or empty selections", {
  expect_error(
    select_query_taxa(df, query_taxa = "not_a_taxon"),
    "not a taxon|not taxa"
  )
  expect_error(
    select_query_taxa(df, query_taxa = ref_taxa, exclude_taxa = ref_taxa),
    "No taxon left"
  )
})

test_that("place_pq accepts a refpkg in place of alignment and tree", {
  skip_if_not_installed("DECIPHER")
  external <- ape::rtree(length(ref_taxa), tip.label = ref_taxa)
  rp <- ref_package_pq(ref_pq, tree = external, model = "GTR+G+I")

  cmds <- place_pq(df, refpkg = rp, exclude_taxa = ref_taxa, cmd_is_run = FALSE)
  expect_named(cmds, c("mafft", "epa-ng"))
  # The model of the package is used unless the caller overrides it
  expect_match(cmds[["epa-ng"]], "GTR\\+G\\+I")
  expect_match(
    place_pq(df, refpkg = rp, model = "JC", cmd_is_run = FALSE)[["epa-ng"]],
    "JC"
  )
})

test_that("place_pq rejects a missing or doubled reference", {
  skip_if_not_installed("DECIPHER")
  expect_error(place_pq(df, cmd_is_run = FALSE), "reference alignment")
  expect_error(place_pq(df, refpkg = "nope", cmd_is_run = FALSE), "refpkg")

  external <- ape::rtree(length(ref_taxa), tip.label = ref_taxa)
  rp <- ref_package_pq(ref_pq, tree = external)
  expect_error(
    place_pq(df, ref_alignment = "a.fasta", refpkg = rp, cmd_is_run = FALSE),
    "not both"
  )
})

test_that("the three reference workflows reach a real placement", {
  skip_if_not(is_epang_installed())
  skip_if_not(MiscMetabar::is_mafft_installed())
  skip_if_not_installed("phangorn")
  skip_if_not_installed("BoSSA")

  # Workflow 2: reference sequences, tree inferred from them
  rp <- ref_package_pq(ref_pq, align_method = "mafft")
  res <- place_pq(df, refpkg = rp, exclude_taxa = ref_taxa)
  expect_setequal(res$query, qry_taxa)
  expect_true(all(res$like_weight_ratio > 0))

  # Workflow 1: the same alignment with an external topology
  rp_ext <- ref_package_pq(
    ref_pq,
    tree = ape::rtree(length(ref_taxa), tip.label = ref_taxa),
    align_method = "mafft"
  )
  res_ext <- place_pq(df, refpkg = rp_ext, query_taxa = qry_taxa)
  expect_setequal(res_ext$query, qry_taxa)

  # Workflow 3: two reference taxa moved into the reference are not placed
  own <- ref_taxa[1:2]
  res3 <- place_pq(df, refpkg = rp, exclude_taxa = c(ref_taxa, qry_taxa[1]))
  expect_false(any(own %in% res3$query))
  expect_false(qry_taxa[1] %in% res3$query)
})
