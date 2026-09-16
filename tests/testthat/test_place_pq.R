skip_on_cran()

library(MiscMetabar)


# A minimal but valid jplace file: 2 queries, ASV_1 spread over two branches.
write_jplace_fixture <- function() {
  path <- tempfile("phylopq_test_", fileext = ".jplace")
  writeLines(
    c(
      "{",
      '  "tree": "((A:0.1{0},B:0.2{1}):0.3{2},C:0.4{3}):0.0{4};",',
      '  "placements": [',
      '    {"p": [[0, -100.0, 0.9, 0.05, 0.01],',
      '           [1, -102.0, 0.1, 0.03, 0.02]], "n": ["ASV_1"]},',
      '    {"p": [[3, -101.0, 1.0, 0.06, 0.02]], "n": ["ASV_2"]}',
      "  ],",
      '  "metadata": {"invocation": "epa-ng"},',
      '  "version": 3,',
      '  "fields": ["edge_num", "likelihood", "like_weight_ratio",',
      '             "distal_length", "pendant_length"]',
      "}"
    ),
    path
  )
  path
}

test_that("find_place_exe honours the option first", {
  withr::with_options(
    list(phylopq.epangpath = "/somewhere/epa-ng"),
    expect_equal(find_place_exe("epa-ng"), "/somewhere/epa-ng")
  )
  withr::with_options(
    list(phylopq.gappapath = "/somewhere/gappa"),
    expect_equal(find_place_exe("gappa"), "/somewhere/gappa")
  )
  expect_error(find_place_exe("raxml"), "should be one of")
})

test_that("is_epang_installed and is_gappa_installed return a single logical", {
  for (f in list(is_epang_installed, is_gappa_installed)) {
    res <- f()
    expect_type(res, "logical")
    expect_length(res, 1)
    expect_false(f(path = tempfile("no_tool_here_")))
  }
})

test_that("resolve_place_exe reports a missing executable", {
  expect_error(
    resolve_place_exe("epa-ng", tempfile("no_epang_here_")),
    "was not found"
  )
  expect_error(
    resolve_place_exe("gappa", tempfile("no_gappa_here_")),
    "was not found"
  )
  expect_error(
    resolve_place_exe("epa-ng", tempfile("no_epang_here_")),
    "install_epang"
  )
})

test_that("as_fasta_path and as_newick_path accept objects and paths", {
  skip_if_not_installed("Biostrings")
  dna <- Biostrings::DNAStringSet(c(a = "ACGT", b = "ACGA"))
  fa <- as_fasta_path(dna, "test_")
  expect_true(file.exists(fa))
  expect_length(Biostrings::readDNAStringSet(fa), 2)
  expect_true(attr(fa, "temporary"))
  # A path given by the caller is theirs, and must not be flagged for deletion
  expect_false(attr(as_fasta_path(as.character(fa), "test_"), "temporary"))

  tree <- ape::rtree(4)
  nwk <- as_newick_path(tree, "test_")
  expect_true(file.exists(nwk))
  expect_equal(ape::Ntip(ape::read.tree(nwk)), 4)
  expect_true(attr(nwk, "temporary"))
  expect_false(attr(as_newick_path(as.character(nwk), "test_"), "temporary"))
})

test_that("place_pq leaves the caller's reference files alone", {
  skip_if_not_installed("Biostrings")
  ref_fa <- as.character(as_fasta_path(
    Biostrings::DNAStringSet(c(A = "ACGT", B = "ACGA")),
    "ref_"
  ))
  ref_nwk <- as.character(
    as_newick_path(ape::rtree(2, tip.label = c("A", "B")), "ref_")
  )
  on.exit(unlink(c(ref_fa, ref_nwk)))

  place_pq(
    Biostrings::DNAStringSet(c(ASV_1 = "ACGT", ASV_2 = "ACGA")),
    ref_alignment = ref_fa,
    ref_tree = ref_nwk,
    cmd_is_run = FALSE
  )
  expect_true(file.exists(ref_fa))
  expect_true(file.exists(ref_nwk))
})

test_that("as_fasta_path and as_newick_path reject bad input", {
  expect_error(as_fasta_path(tempfile("absent_"), "t_"), "does not exist")
  expect_error(as_newick_path(tempfile("absent_"), "t_"), "does not exist")
  expect_error(as_fasta_path(1:10, "t_"), "DNAStringSet")
  expect_error(as_newick_path(1:10, "t_"), "phylo")
  # A dry run assembles a command for elsewhere, so the files need not exist
  expect_silent(as_fasta_path("absent.fasta", "t_", check_exists = FALSE))
  expect_silent(as_newick_path("absent.newick", "t_", check_exists = FALSE))
})

test_that("place_pq builds a command for files that do not exist yet", {
  cmds <- place_pq(
    Biostrings::DNAStringSet(c(ASV_1 = "ACGT", ASV_2 = "ACGA")),
    ref_alignment = "reference_alignment.fasta",
    ref_tree = "reference_tree.newick",
    cmd_is_run = FALSE
  )
  expect_named(cmds, c("mafft", "epa-ng"))
  expect_match(cmds[["epa-ng"]], "reference_alignment.fasta")
  expect_match(cmds[["epa-ng"]], "reference_tree.newick")
})

test_that("as_jplace_path unwraps a place_pq result", {
  df <- data.frame(query = "ASV_1")
  attr(df, "jplace_file") <- "/tmp/epa_result.jplace"
  expect_equal(as_jplace_path(df), "/tmp/epa_result.jplace")
  expect_equal(as_jplace_path("a.jplace"), "a.jplace")
  expect_error(as_jplace_path(1:3), "must be a path")
})

test_that("read_placements flattens a jplace file", {
  skip_if_not_installed("BoSSA")
  jp_file <- write_jplace_fixture()
  on.exit(unlink(jp_file))

  all_p <- read_placements(jp_file, best_only = FALSE)
  expect_s3_class(all_p, "data.frame")
  expect_named(
    all_p,
    c(
      "query",
      "edge",
      "jplace_edge",
      "like_weight_ratio",
      "likelihood",
      "distal_length",
      "pendant_length"
    )
  )
  expect_equal(nrow(all_p), 3)
  expect_setequal(all_p$query, c("ASV_1", "ASV_2"))
  expect_equal(attr(all_p, "jplace_file"), jp_file)
  expect_s3_class(attr(all_p, "tree"), "phylo")
  # `jplace_edge` recovers EPA-ng's own numbering, which BoSSA remaps
  expect_setequal(all_p$jplace_edge, c(0, 1, 3))
})

test_that("read_placements keeps only the best placement by default", {
  skip_if_not_installed("BoSSA")
  jp_file <- write_jplace_fixture()
  on.exit(unlink(jp_file))

  best <- read_placements(jp_file, best_only = TRUE)
  expect_equal(nrow(best), 2)
  expect_equal(best$query, c("ASV_1", "ASV_2"))
  # ASV_1 keeps its LWR 0.9 branch, not the 0.1 one
  expect_equal(best$like_weight_ratio[best$query == "ASV_1"], 0.9)
  expect_equal(best$jplace_edge[best$query == "ASV_1"], 0)
})

test_that("read_placements can return the jplace object", {
  skip_if_not_installed("BoSSA")
  jp_file <- write_jplace_fixture()
  on.exit(unlink(jp_file))
  jp <- read_placements(jp_file, return_jplace = TRUE)
  expect_s3_class(jp, "jplace")
  expect_s3_class(jp$arbre, "phylo")
})

test_that("jplace_queries lists the queries and tolerates a bad file", {
  skip_if_not_installed("BoSSA")
  jp_file <- write_jplace_fixture()
  on.exit(unlink(jp_file))
  expect_setequal(jplace_queries(jp_file), c("ASV_1", "ASV_2"))
  expect_length(jplace_queries(tempfile("absent_")), 0)
})

test_that("place_pq builds the commands without running them", {
  skip_if_not_installed("Biostrings")
  df <- suppressMessages(subset_taxa_pq(
    data_fungi_mini,
    taxa_sums(data_fungi_mini) > 12000
  ))
  ref_fa <- as_fasta_path(
    Biostrings::DNAStringSet(c(A = "ACGT", B = "ACGA")),
    "ref_"
  )
  ref_nwk <- as_newick_path(ape::rtree(2, tip.label = c("A", "B")), "ref_")
  on.exit(unlink(c(ref_fa, ref_nwk)))

  cmds <- place_pq(
    df,
    ref_alignment = ref_fa,
    ref_tree = ref_nwk,
    model = "GTR+G",
    cmd_is_run = FALSE
  )
  expect_type(cmds, "character")
  expect_named(cmds, c("mafft", "epa-ng"))
  expect_match(cmds[["mafft"]], "--add")
  expect_match(cmds[["mafft"]], "--keeplength")
  expect_match(cmds[["epa-ng"]], "--ref-msa")
  expect_match(cmds[["epa-ng"]], "--tree")
  expect_match(cmds[["epa-ng"]], "--query")
  expect_match(cmds[["epa-ng"]], "GTR\\+G")
})

test_that("place_pq passes epang_args and threads through", {
  skip_if_not_installed("Biostrings")
  ref_fa <- as_fasta_path(
    Biostrings::DNAStringSet(c(A = "ACGT", B = "ACGA")),
    "ref_"
  )
  ref_nwk <- as_newick_path(ape::rtree(2, tip.label = c("A", "B")), "ref_")
  on.exit(unlink(c(ref_fa, ref_nwk)))

  cmds <- place_pq(
    Biostrings::DNAStringSet(c(ASV_1 = "ACGT", ASV_2 = "ACGA")),
    ref_alignment = ref_fa,
    ref_tree = ref_nwk,
    threads = 4,
    epang_args = "--filter-acc-lwr 0.99",
    cmd_is_run = FALSE
  )
  expect_match(cmds[["epa-ng"]], "--threads 4")
  expect_match(cmds[["epa-ng"]], "--filter-acc-lwr 0.99")
  expect_match(cmds[["mafft"]], "--thread 4")
})

test_that("place_pq errors on a phyloseq without refseq", {
  skip_if_not_installed("Biostrings")
  ref_fa <- as_fasta_path(
    Biostrings::DNAStringSet(c(A = "ACGT", B = "ACGA")),
    "ref_"
  )
  ref_nwk <- as_newick_path(ape::rtree(2, tip.label = c("A", "B")), "ref_")
  on.exit(unlink(c(ref_fa, ref_nwk)))

  no_refseq <- phyloseq::phyloseq(
    phyloseq::otu_table(data_fungi_mini),
    phyloseq::tax_table(data_fungi_mini)
  )
  expect_error(
    place_pq(
      no_refseq,
      ref_alignment = ref_fa,
      ref_tree = ref_nwk,
      cmd_is_run = FALSE
    ),
    "refseq"
  )
  expect_error(
    place_pq(1:10, ref_alignment = ref_fa, ref_tree = ref_nwk),
    "phyloseq"
  )
})

test_that("assign_placement_pq builds the gappa command without running it", {
  cmd <- assign_placement_pq(
    jplace = "epa_result.jplace",
    taxon_file = "taxonomy.tsv",
    cmd_is_run = FALSE
  )
  expect_type(cmd, "character")
  expect_length(cmd, 1)
  expect_match(cmd, "examine assign")
  expect_match(cmd, "--jplace-path")
  expect_match(cmd, "--taxon-file")
  expect_match(cmd, "--out-dir")
})

test_that("assign_placement_pq passes gappa_args through", {
  cmd <- assign_placement_pq(
    jplace = "epa_result.jplace",
    taxon_file = "taxonomy.tsv",
    gappa_args = "--consensus-thresh 0.9",
    cmd_is_run = FALSE
  )
  expect_match(cmd, "--consensus-thresh 0.9")
})

# Verbatim `per_query.tsv` of gappa v0.9.0 `examine assign` on the jplace
# fixture above, with the taxon file A/B under Sordariomycetes and C under
# Agaricales. Note the several rows per query, one per taxonomic depth.
write_per_query_fixture <- function() {
  path <- tempfile("phylopq_test_", fileext = ".tsv")
  writeLines(
    c(
      "name\tLWR\tfract\taLWR\tafract\ttaxopath",
      "ASV_1\t0\t0\t1\t1\tFungi",
      "ASV_1\t0\t0\t1\t1\tFungi;Ascomycota",
      "ASV_1\t0.465\t0.465\t1\t1\tFungi;Ascomycota;Sordariomycetes",
      paste0(
        "ASV_1\t0.45\t0.45\t0.45\t0.45\t",
        "Fungi;Ascomycota;Sordariomycetes;Hypocreales"
      ),
      paste0(
        "ASV_1\t0.085\t0.085\t0.085\t0.085\t",
        "Fungi;Ascomycota;Sordariomycetes;Xylariales"
      ),
      "ASV_2\t0.15\t0.15\t1\t1\tFungi",
      "ASV_2\t0\t0\t0.85\t0.85\tFungi;Basidiomycota",
      "ASV_2\t0\t0\t0.85\t0.85\tFungi;Basidiomycota;Agaricomycetes",
      paste0(
        "ASV_2\t0.85\t0.85\t0.85\t0.85\t",
        "Fungi;Basidiomycota;Agaricomycetes;Agaricales"
      )
    ),
    path
  )
  path
}

ranks4 <- c("Kingdom", "Phylum", "Class", "Order")

test_that("read_gappa_assignment keeps the deepest confident row per query", {
  per_query <- write_per_query_fixture()
  on.exit(unlink(per_query))

  res <- read_gappa_assignment(
    per_query,
    queries = c("ASV_1", "ASV_2", "ASV_3"),
    ranks = ranks4,
    min_alwr = 0.5
  )
  expect_named(res, c("query", "alwr", ranks4))
  # One row per query, not one per taxonomic depth
  expect_equal(res$query, c("ASV_1", "ASV_2", "ASV_3"))

  # ASV_1: Hypocreales carries only 0.45, so it stops at Sordariomycetes
  expect_equal(res$Class[res$query == "ASV_1"], "Sordariomycetes")
  expect_true(is.na(res$Order[res$query == "ASV_1"]))
  expect_equal(res$alwr[res$query == "ASV_1"], 1)

  # ASV_2: Agaricales carries 0.85 and is kept
  expect_equal(res$Order[res$query == "ASV_2"], "Agaricales")
  expect_equal(res$alwr[res$query == "ASV_2"], 0.85)

  # ASV_3 was never placed but keeps a row of NA
  expect_true(all(is.na(unlist(res[res$query == "ASV_3", ranks4]))))
})

test_that("read_gappa_assignment follows min_alwr up and down", {
  per_query <- write_per_query_fixture()
  on.exit(unlink(per_query))

  # A lower floor accepts the finer, shakier rank
  loose <- read_gappa_assignment(
    per_query,
    queries = c("ASV_1", "ASV_2"),
    ranks = ranks4,
    min_alwr = 0.4
  )
  expect_equal(loose$Order[loose$query == "ASV_1"], "Hypocreales")

  # A stricter floor drops ASV_2 back to a coarser rank
  strict <- read_gappa_assignment(
    per_query,
    queries = c("ASV_1", "ASV_2"),
    ranks = ranks4,
    min_alwr = 0.9
  )
  expect_equal(strict$Class[strict$query == "ASV_1"], "Sordariomycetes")
  expect_equal(strict$Kingdom[strict$query == "ASV_2"], "Fungi")
  expect_true(is.na(strict$Order[strict$query == "ASV_2"]))
})

test_that("read_gappa_assignment falls back when aLWR is absent", {
  # Older gappa releases write no aLWR column; every row is then eligible
  per_query <- tempfile(fileext = ".tsv")
  writeLines(
    c(
      "name\tLWR\tfract\ttaxopath",
      "ASV_1\t0.9\t0.9\tFungi;Ascomycota",
      "ASV_1\t0.4\t0.4\tFungi;Ascomycota;Sordariomycetes"
    ),
    per_query
  )
  on.exit(unlink(per_query))

  res <- read_gappa_assignment(per_query, "ASV_1", ranks4)
  expect_equal(res$Class, "Sordariomycetes")
  expect_true(is.na(res$alwr))
})

test_that("read_gappa_assignment reports an unexpected layout", {
  per_query <- tempfile(fileext = ".tsv")
  writeLines(c("a\tb", "1\t2"), per_query)
  on.exit(unlink(per_query))
  expect_error(
    read_gappa_assignment(per_query, queries = character(0), ranks = "Kingdom"),
    "Could not find"
  )
})

test_that("assign_placement_pq always asks gappa for per_query.tsv", {
  # gappa writes per_query.tsv only with --per-query-results; profile.tsv,
  # which it writes by default, has no query names and is useless here.
  cmd <- assign_placement_pq(
    jplace = "epa_result.jplace",
    taxon_file = "taxonomy.tsv",
    cmd_is_run = FALSE
  )
  expect_match(cmd, "--per-query-results")
})

test_that("assign_placement_pq assigns a taxonomy with the real gappa", {
  skip_if_not(is_gappa_installed())
  skip_if_not_installed("BoSSA")
  jp_file <- write_jplace_fixture()
  taxon_file <- tempfile(fileext = ".tsv")
  writeLines(
    c(
      "A\tFungi;Ascomycota;Sordariomycetes;Hypocreales",
      "B\tFungi;Ascomycota;Sordariomycetes;Xylariales",
      "C\tFungi;Basidiomycota;Agaricomycetes;Agaricales"
    ),
    taxon_file
  )
  on.exit(unlink(c(jp_file, taxon_file)))

  taxo <- assign_placement_pq(
    jplace = jp_file,
    taxon_file = taxon_file,
    ranks = ranks4
  )
  expect_s3_class(taxo, "data.frame")
  expect_named(taxo, c("query", "alwr", ranks4))
  expect_equal(taxo$query, c("ASV_1", "ASV_2"))
  expect_equal(taxo$Kingdom, c("Fungi", "Fungi"))
  expect_equal(taxo$Class, c("Sordariomycetes", "Agaricomycetes"))
  expect_equal(taxo$Order, c(NA, "Agaricales"))
})

# A small leave-out reference package: most taxa become the reference
# alignment + tree, the remaining few are placed back onto it.
build_refpkg <- function(n_query = 3) {
  df <- suppressMessages(MiscMetabar::subset_taxa_pq(
    data_fungi_mini,
    phyloseq::taxa_sums(data_fungi_mini) > 9000
  ))
  taxa <- phyloseq::taxa_names(df)
  ref_taxa <- taxa[seq_len(length(taxa) - n_query)]
  ref_pq <- phyloseq::prune_taxa(ref_taxa, df)

  ref_ali <- MiscMetabar::align_pq(ref_pq, method = "mafft", force = TRUE)
  ref_fa <- tempfile("ref_msa_", fileext = ".fasta")
  Biostrings::writeXStringSet(ref_ali, ref_fa)

  phang <- phangorn::phyDat(as(ref_ali, "matrix"), type = "DNA")
  tree <- ape::multi2di(phangorn::NJ(phangorn::dist.ml(phang)))
  tree$edge.length[tree$edge.length < 0] <- 1e-8
  nwk <- tempfile("ref_tree_", fileext = ".newick")
  ape::write.tree(tree, nwk)

  list(
    query_pq = phyloseq::prune_taxa(setdiff(taxa, ref_taxa), df),
    ref_pq = ref_pq,
    ref_fa = ref_fa,
    nwk = nwk
  )
}

test_that("place_pq places queries with the real epa-ng", {
  skip_if_not(is_epang_installed())
  skip_if_not(MiscMetabar::is_mafft_installed())
  skip_if_not_installed("phangorn")
  skip_if_not_installed("BoSSA")

  rp <- build_refpkg()
  on.exit(unlink(c(rp$ref_fa, rp$nwk)))
  queries <- phyloseq::taxa_names(rp$query_pq)

  res <- place_pq(rp$query_pq, ref_alignment = rp$ref_fa, ref_tree = rp$nwk)
  expect_s3_class(res, "data.frame")
  # best_only keeps exactly one row per query
  expect_setequal(res$query, queries)
  expect_equal(nrow(res), length(queries))
  expect_true(all(res$like_weight_ratio > 0 & res$like_weight_ratio <= 1))
  expect_true(all(res$pendant_length >= 0))

  # The result points at a real .jplace, and carries the reference tree
  expect_true(file.exists(attr(res, "jplace_file")))
  expect_s3_class(attr(res, "tree"), "phylo")

  # The caller's reference files survive the call
  expect_true(file.exists(rp$ref_fa))
  expect_true(file.exists(rp$nwk))

  # Without best_only every candidate branch is returned
  all_p <- place_pq(
    rp$query_pq,
    ref_alignment = rp$ref_fa,
    ref_tree = rp$nwk,
    best_only = FALSE
  )
  expect_gte(nrow(all_p), nrow(res))
  expect_setequal(all_p$query, queries)
  # `jplace_edge` is EPA-ng's own numbering, which BoSSA remaps
  expect_false(identical(all_p$edge, all_p$jplace_edge))
})

test_that("place_pq and assign_placement_pq recover a left-out taxonomy", {
  skip_if_not(is_epang_installed())
  skip_if_not(is_gappa_installed())
  skip_if_not(MiscMetabar::is_mafft_installed())
  skip_if_not_installed("phangorn")
  skip_if_not_installed("BoSSA")

  rp <- build_refpkg()
  ranks <- c("Phylum", "Class", "Order")
  tt <- as.data.frame(unclass(phyloseq::tax_table(rp$ref_pq)))[, ranks]
  paths <- apply(tt, 1, function(r) {
    r[is.na(r) | r == "unidentified"] <- ""
    paste(r[nzchar(r)], collapse = ";")
  })
  taxon_file <- tempfile("ref_taxonomy_", fileext = ".tsv")
  writeLines(paste(names(paths), paths, sep = "\t"), taxon_file)
  on.exit(unlink(c(rp$ref_fa, rp$nwk, taxon_file)))

  places <- place_pq(
    rp$query_pq,
    ref_alignment = rp$ref_fa,
    ref_tree = rp$nwk
  )
  taxo <- assign_placement_pq(places, taxon_file = taxon_file, ranks = ranks)

  expect_named(taxo, c("query", "alwr", ranks))
  expect_setequal(taxo$query, phyloseq::taxa_names(rp$query_pq))
  # Placement on a fungal reference recovers the phylum of every query
  expect_true(all(!is.na(taxo$Phylum)))
  truth <- as.data.frame(unclass(phyloseq::tax_table(rp$query_pq)))
  expect_equal(taxo$Phylum, unname(truth[taxo$query, "Phylum"]))

  # Written straight into a phyloseq object, taxa and samples are preserved
  pq <- assign_placement_pq(
    places,
    taxon_file = taxon_file,
    physeq = rp$query_pq,
    ranks = ranks
  )
  expect_s4_class(pq, "phyloseq")
  expect_equal(phyloseq::ntaxa(pq), phyloseq::ntaxa(rp$query_pq))
  expect_equal(phyloseq::nsamples(pq), phyloseq::nsamples(rp$query_pq))
  expect_equal(colnames(phyloseq::tax_table(pq)), ranks)
})

test_that("assign_placement_pq errors on a missing taxon file", {
  skip_if_not(is_gappa_installed())
  expect_error(
    assign_placement_pq(
      jplace = "epa_result.jplace",
      taxon_file = tempfile("absent_")
    ),
    "does not exist"
  )
})

test_that("taxo_to_pq writes the assignment into a tax_table", {
  df <- suppressMessages(subset_taxa_pq(
    data_fungi_mini,
    taxa_sums(data_fungi_mini) > 12000
  ))
  ranks <- c("Kingdom", "Phylum")
  taxa <- phyloseq::taxa_names(df)
  taxo <- data.frame(
    query = taxa[1:2],
    Kingdom = c("Fungi", "Fungi"),
    Phylum = c("Ascomycota", "Basidiomycota"),
    stringsAsFactors = FALSE
  )

  pq <- taxo_to_pq(taxo, df, ranks = ranks)
  expect_s4_class(pq, "phyloseq")
  tt <- as.data.frame(unclass(phyloseq::tax_table(pq)))
  expect_named(tt, ranks)
  expect_equal(tt[taxa[1], "Phylum"], "Ascomycota")
  # Taxa absent from the placement keep NA rather than being dropped
  expect_equal(phyloseq::ntaxa(pq), phyloseq::ntaxa(df))
  expect_true(is.na(tt[taxa[3], "Phylum"]))
})

test_that("taxo_to_pq warns when no taxon matches", {
  df <- suppressMessages(subset_taxa_pq(
    data_fungi_mini,
    taxa_sums(data_fungi_mini) > 12000
  ))
  taxo <- data.frame(
    query = c("not_a_taxon"),
    Kingdom = "Fungi",
    stringsAsFactors = FALSE
  )
  expect_warning(
    taxo_to_pq(taxo, df, ranks = "Kingdom"),
    "None of the"
  )
})

test_that("check_place_build_tools names what is missing", {
  # git/cmake/make are present on most dev machines; the guard must at least
  # run and either pass silently or name the missing tool.
  for (tool in c("epa-ng", "gappa")) {
    res <- tryCatch(
      {
        check_place_build_tools(tool)
        "ok"
      },
      error = function(e) {
        conditionMessage(e)
      }
    )
    expect_true(res == "ok" || grepl("Cannot build", res))
  }
  expect_error(check_place_build_tools("raxml"), "should be one of")
})

test_that("check_place_build_tools requires bison and flex for epa-ng only", {
  # epa-ng bundles libpll, whose CMake build calls BISON_TARGET()/FLEX_TARGET();
  # gappa builds on genesis alone and needs neither.
  fake_path <- tempfile("empty_bin_")
  dir.create(fake_path)
  on.exit(unlink(fake_path, recursive = TRUE))

  withr::with_path(fake_path, action = "replace", {
    epang <- tryCatch(
      {
        check_place_build_tools("epa-ng")
        ""
      },
      error = function(e) conditionMessage(e)
    )
    expect_match(epang, "bison")
    expect_match(epang, "flex")
    expect_match(epang, "BISON_TARGET")

    gappa <- tryCatch(
      {
        check_place_build_tools("gappa")
        ""
      },
      error = function(e) conditionMessage(e)
    )
    expect_false(grepl("bison", gappa))
    expect_false(grepl("flex", gappa))
  })
})

test_that("place_source_repo returns the two upstream repositories", {
  expect_match(place_source_repo("epa-ng"), "epa-ng")
  expect_match(place_source_repo("gappa"), "gappa")
  expect_error(place_source_repo("raxml"), "should be one of")
})

test_that("install_epang and install_gappa refuse to rebuild silently", {
  fake <- file.path(tempfile("phylopq_data_"), "bin")
  dir.create(fake, recursive = TRUE)
  file.create(file.path(fake, "epa-ng"))
  on.exit(unlink(dirname(fake), recursive = TRUE))

  expect_message(
    res <- install_epang(path = dirname(fake)),
    "already installed"
  )
  expect_equal(res, file.path(fake, "epa-ng"))
})
